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
/// The material LUT lives in `aod::GeomCCDBObjects` rather than here: it belongs to
/// the geometry/material family, whose validity is essentially static, and keeping it
/// out means joining `aod::GloCCDBObjects` does not drag in a multi-hundred-MB object
/// nobody asked for. Join whichever tables you need — the duplicated `aod::Timestamps`
/// is deduplicated:
/// \code
///   using BCsWithLUT = soa::Join<aod::BCsWithTimestamps, aod::GeomCCDBObjects>;
///   // rectifyPtrFromFile() is applied by the column's finaliser, so the
///   // object handed back is ready to use; the task only has to (re)install it.
///   auto* lut = &bcs.begin().matLUT();
///   if (lut != mLastLUT) {
///     o2::base::Propagator::Instance()->setMatLUT(lut);
///     mLastLUT = lut;
///   }
/// \endcode

#ifndef COMMON_DATAMODEL_GLOCCDBOBJECTS_H_
#define COMMON_DATAMODEL_GLOCCDBOBJECTS_H_

#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPECSObject.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
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

/// The material LUT is a FlatObject: straight out of the ROOT streamer its internal
/// pointers are unfixed and its voxel lookup is unbuilt, so it is finalised with
/// MatLayerCylSet::rectifyPtrFromFile() before ever being handed to a task.
DECLARE_SOA_CCDB_COLUMN_FULL(MatLUT, "fMatLUT", matLUT, o2::base::MatLayerCylSet, "GLO/Param/MatLUT", //!
                             [](o2::base::MatLayerCylSet* lut) { return o2::base::MatLayerCylSet::rectifyPtrFromFile(lut); });
} // namespace ccdbGlo

/// Full table — join with aod::BCsWithTimestamps to obtain all four objects.
DECLARE_SOA_TIMESTAMPED_TABLE(GloCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "GLOCCDBOBJ", //!
                              ccdbGlo::GRPMagField, ccdbGlo::MeanVertex, ccdbGlo::GRPECSObject, ccdbGlo::GRPLHCIFData);

/// Geometry and material description: objects which describe where the detector material
/// is, and which share an essentially static interval of validity. Kept apart from the GRP
/// family above, which changes per run (and, for GRPMagField, per timeframe).
/// The aligned/ideal geometry and the per-detector alignment objects belong here too when
/// they get columns; see GRPGeomRequest in O2 (GLO/Config/GeometryAligned, GLO/Config/Geometry,
/// <DET>/Calib/Align) for the family.
/// Join it alongside aod::GloCCDBObjects when a task needs both — the duplicated
/// aod::Timestamps is deduplicated when the joined table's originals are merged.
/// Uniform in the run number: the geometry/material description does not change within a
/// run, so the fetcher queries once per distinct run rather than once per BC.
DECLARE_SOA_UNIFORM_TABLE(GeomCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp,
                          aod::BCs, o2::aod::bc::RunNumber, 1, "GEOMCCDBOBJ", //!
                          ccdbGlo::MatLUT);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_GLOCCDBOBJECTS_H_
