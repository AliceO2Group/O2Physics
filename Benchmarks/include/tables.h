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

DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax0R0, realQMax0R0, float, "fRealQMax0R0");
DECLARE_SOA_COLUMN(RealQMax0R0C, realQMax0R0C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax1R0, realQMax1R0, float, "fRealQMax1R0");
DECLARE_SOA_COLUMN(RealQMax1R0C, realQMax1R0C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax2R0, realQMax2R0, float, "fRealQMax2R0");
DECLARE_SOA_COLUMN(RealQMax2R0C, realQMax2R0C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax3R0, realQMax3R0, float, "fRealQMax3R0");
DECLARE_SOA_COLUMN(RealQMax3R0C, realQMax3R0C, float);

DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot0R0, realQTot0R0, float, "fRealQTot0R0");
DECLARE_SOA_COLUMN(RealQTot0R0C, realQTot0R0C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot1R0, realQTot1R0, float, "fRealQTot1R0");
DECLARE_SOA_COLUMN(RealQTot1R0C, realQTot1R0C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot2R0, realQTot2R0, float, "fRealQTot2R0");
DECLARE_SOA_COLUMN(RealQTot2R0C, realQTot2R0C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot3R0, realQTot3R0, float, "fRealQTot3R0");
DECLARE_SOA_COLUMN(RealQTot3R0C, realQTot3R0C, float);

DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot0, realQMaxTot0, float, "fRealQMaxTot0");
DECLARE_SOA_COLUMN(RealQMaxTot0C, realQMaxTot0C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot1, realQMaxTot1, float, "fRealQMaxTot1");
DECLARE_SOA_COLUMN(RealQMaxTot1C, realQMaxTot1C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot2, realQMaxTot2, float, "fRealQMaxTot2");
DECLARE_SOA_COLUMN(RealQMaxTot2C, realQMaxTot2C, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot3, realQMaxTot3, float, "fRealQMaxTot3");
DECLARE_SOA_COLUMN(RealQMaxTot3C, realQMaxTot3C, float);

DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax0R0_mad, realQMax0R0_mad, float, "fRealQMax0R0_mad");
DECLARE_SOA_COLUMN(RealQMax0R0C_mad, realQMax0R0C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax1R0_mad, realQMax1R0_mad, float, "fRealQMax1R0_mad");
DECLARE_SOA_COLUMN(RealQMax1R0C_mad, realQMax1R0C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax2R0_mad, realQMax2R0_mad, float, "fRealQMax2R0_mad");
DECLARE_SOA_COLUMN(RealQMax2R0C_mad, realQMax2R0C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMax3R0_mad, realQMax3R0_mad, float, "fRealQMax3R0_mad");
DECLARE_SOA_COLUMN(RealQMax3R0C_mad, realQMax3R0C_mad, float);

DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot0R0_mad, realQTot0R0_mad, float, "fRealQTot0R0_mad");
DECLARE_SOA_COLUMN(RealQTot0R0C_mad, realQTot0R0C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot1R0_mad, realQTot1R0_mad, float, "fRealQTot1R0_mad");
DECLARE_SOA_COLUMN(RealQTot1R0C_mad, realQTot1R0C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot2R0_mad, realQTot2R0_mad, float, "fRealQTot2R0_mad");
DECLARE_SOA_COLUMN(RealQTot2R0C_mad, realQTot2R0C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQTot3R0_mad, realQTot3R0_mad, float, "fRealQTot3R0_mad");
DECLARE_SOA_COLUMN(RealQTot3R0C_mad, realQTot3R0C_mad, float);

DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot0_mad, realQMaxTot0_mad, float, "fRealQMaxTot0_mad");
DECLARE_SOA_COLUMN(RealQMaxTot0C_mad, realQMaxTot0C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot1_mad, realQMaxTot1_mad, float, "fRealQMaxTot1_mad");
DECLARE_SOA_COLUMN(RealQMaxTot1C_mad, realQMaxTot1C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot2_mad, realQMaxTot2_mad, float, "fRealQMaxTot2_mad");
DECLARE_SOA_COLUMN(RealQMaxTot2C_mad, realQMaxTot2C_mad, float);
DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(RealQMaxTot3_mad, realQMaxTot3_mad, float, "fRealQMaxTot3_mad");
DECLARE_SOA_COLUMN(RealQMaxTot3C_mad, realQMaxTot3C_mad, float);

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
DECLARE_SOA_COLUMN(Correction1, correction1, float);
}

// intermediate table
DECLARE_SOA_TABLE(TracksTemporaryExtra, "AOD", "TRKTEMPEX",
                  intermediate::HRate, intermediate::ClampedTPCMult, intermediate::Occupancy, intermediate::Correction1,
                  aod::track::TPCSignal, aod::track::Signed1Pt, aod::track::Tgl);
using TracksQAEx = soa::Join<TracksQAVersion, TracksTemporaryExtra>;

// final table
DECLARE_SOA_CONFIGURABLE_EXTENDED_TABLE(TracksQACorrectedE, TracksQAEx, "TRKQACORE", extensions::RealTPCSignalN);
using MoreTracksFinal = soa::Join<TracksQAVersion, TracksQACorrectedECfgExtension>;

DECLARE_SOA_CONFIGURABLE_EXTENDED_TABLE(TracksQACorrectedEFull, TracksQAEx, "TRKQACOREF",
                                        extensions::RealTPCSignalN,
                                        extensions::RealQMax0R0, extensions::RealQMax1R0, extensions::RealQMax2R0, extensions::RealQMax3R0,
                                        extensions::RealQMax0R0_mad, extensions::RealQMax1R0_mad, extensions::RealQMax2R0_mad, extensions::RealQMax3R0_mad,
                                        extensions::RealQTot0R0, extensions::RealQTot1R0, extensions::RealQTot2R0, extensions::RealQTot3R0,
                                        extensions::RealQTot0R0_mad, extensions::RealQTot1R0_mad, extensions::RealQTot2R0_mad, extensions::RealQTot3R0_mad,
                                        extensions::RealQMaxTot0, extensions::RealQMaxTot1, extensions::RealQMaxTot2, extensions::RealQMaxTot3,
                                        extensions::RealQMaxTot0_mad, extensions::RealQMaxTot1_mad, extensions::RealQMaxTot2_mad, extensions::RealQMaxTot3_mad);
using MoreTracksFinalFull = soa::Join<TracksQAVersion, TracksQACorrectedEFullCfgExtension>;

// final table for direct calculation
DECLARE_SOA_TABLE(TracksQACorrected, "AOD", "TRKQACOR", extensions::RealTPCSignalNC);
DECLARE_SOA_TABLE(TracksQACorrectedFull, "AOD", "TRKQACORF",
                  extensions::RealTPCSignalNC,
                  extensions::RealQMax0R0C, extensions::RealQMax1R0C, extensions::RealQMax2R0C, extensions::RealQMax3R0C,
                  extensions::RealQMax0R0C_mad, extensions::RealQMax1R0C_mad, extensions::RealQMax2R0C_mad, extensions::RealQMax3R0C_mad,
                  extensions::RealQTot0R0C, extensions::RealQTot1R0C, extensions::RealQTot2R0C, extensions::RealQTot3R0C,
                  extensions::RealQTot0R0C_mad, extensions::RealQTot1R0C_mad, extensions::RealQTot2R0C_mad, extensions::RealQTot3R0C_mad,
                  extensions::RealQMaxTot0C, extensions::RealQMaxTot1C, extensions::RealQMaxTot2C, extensions::RealQMaxTot3C,
                  extensions::RealQMaxTot0C_mad, extensions::RealQMaxTot1C_mad, extensions::RealQMaxTot2C_mad, extensions::RealQMaxTot3C_mad);

using TracksD = soa::Join<TracksIU, ExtTracksD>;
using TracksID = soa::Join<TracksIU, ExtTracksID>;
} // namespace o2::aod
#endif
