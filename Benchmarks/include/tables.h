// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef O2PHYSICS_BENCHMARKS_INCLUDE_TABLES_H_
#define O2PHYSICS_BENCHMARKS_INCLUDE_TABLES_H_
#include <Framework/AnalysisDataModel.h>
#include <Common/DataModel/TrackSelectionTables.h>

namespace o2::aod
{
namespace extensions
{
// example generic functions
// 1. simple dynamic column and a column to be filled with the result
DECLARE_SOA_DYNAMIC_COLUMN(Direct, direct, [](float x, float y, float z, float t) -> float { return t * (x * x + y * y + z * z); });
DECLARE_SOA_COLUMN(DirectM, directm, float);

// 2. arbitrary function dynamic column and a column to be filled with the result
DECLARE_SOA_DYNAMIC_COLUMN(Indirect, indirect, [](float phi, float x, float y, float z, std::function<float(float, float, float, float)> const& f) -> float { return f(phi, x, y, z); });
DECLARE_SOA_COLUMN(IndirectM, indirectm, float);

// 3. arbitrary expression column placeholder
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(Expr, expr, float, "Expr");

// 4. example realistic values - corrected dE/dx
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealTPCSignalN, realTPCSignalN, float, "fRealTPCSignalN");
DECLARE_SOA_COLUMN(RealTPCSignalNC, realTPCSignalNC, float);
} // namespace extensions

// tables with simple and arbitrary function dynamic columns
DECLARE_SOA_TABLE(ExtTracksD, "AOD", "TRKD", extensions::Direct<aod::track::X, aod::track::Y, aod::track::Z>);
DECLARE_SOA_TABLE(ExtTracksID, "AOD", "TRKID", extensions::Indirect<aod::track::Phi, aod::track::X, aod::track::Y, aod::track::Z>);

// tables to be filled with the results of the above
DECLARE_SOA_TABLE(ExtTracksDM, "AOD", "TRKDM", extensions::DirectM);
DECLARE_SOA_TABLE(ExtTracksIDM, "AOD", "TRKIDM", extensions::IndirectM);

// extended table with an arbitrary expression column
DECLARE_SOA_CONFIGURABLE_EXTENDED_TABLE(TracksE, TracksIU, "TRKE", extensions::Expr);

// intermediate values for the realistic calculation
namespace intermediate {
DECLARE_SOA_COLUMN(HRate, hRate, float);
DECLARE_SOA_COLUMN(ClampedTPCMult, clampedTPCmult, float);
DECLARE_SOA_COLUMN(Occupancy, occupancy, float);
DECLARE_SOA_COLUMN(Correction0, correction0, float);
}

// intermediate table
DECLARE_SOA_TABLE(TracksTemporaryExtra, "AOD", "TRKTEMPEX",
                  intermediate::HRate, intermediate::ClampedTPCMult, intermediate::Occupancy, intermediate::Correction0,
                  aod::track::TPCSignal, aod::track::Signed1Pt, aod::track::Tgl);
using TracksQAEx = soa::Join<TracksQAVersion, TracksTemporaryExtra>;

// final table
DECLARE_SOA_CONFIGURABLE_EXTENDED_TABLE(TracksQACorrectedE, TracksQAEx, "TRKWTPCE", extensions::RealTPCSignalN);
using MoreTracksFinal = soa::Join<TracksQAVersion, TracksQACorrectedECfgExtension>;

// final table for direct calculation
DECLARE_SOA_TABLE(TracksQACorrected, "AOD", "TPCEXT", extensions::RealTPCSignalNC);

using TracksD = soa::Join<TracksIU, ExtTracksD>;
using TracksID = soa::Join<TracksIU, ExtTracksID>;
} // namespace o2::aod
#endif
