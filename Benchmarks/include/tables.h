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
#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace extensions
{
DECLARE_SOA_DYNAMIC_COLUMN(Direct, direct, [](float x, float y, float z, float t) -> float { return t * (x * x + y * y + z * z); });
DECLARE_SOA_COLUMN(DirectM, directm, float);

DECLARE_SOA_DYNAMIC_COLUMN(Indirect, indirect, [](float phi, float x, float y, float z, std::function<float(float, float, float, float)> const& f) -> float { return f(phi, x, y, z); });
DECLARE_SOA_COLUMN(IndirectM, indirectm, float);

DECLARE_SOA_CONFIGURABLE_EXPRESSION_COLUMN(Expr, expr, float, "Expr");
} // namespace extensions

DECLARE_SOA_TABLE(ExtTracksD, "AOD", "TRKD", extensions::Direct<aod::track::X, aod::track::Y, aod::track::Z>);
DECLARE_SOA_TABLE(ExtTracksID, "AOD", "TRKID", extensions::Indirect<aod::track::Phi, aod::track::X, aod::track::Y, aod::track::Z>);

DECLARE_SOA_TABLE(ExtTracksDM, "AOD", "TRKDM", extensions::DirectM);
DECLARE_SOA_TABLE(ExtTracksIDM, "AOD", "TRKIDM", extensions::IndirectM);

DECLARE_SOA_CONFIGURABLE_EXTENDED_TABLE(TracksE, TracksIU, "TRKE", extensions::Expr);

using TracksD = soa::Join<TracksIU, ExtTracksD>;
using TracksID = soa::Join<TracksIU, ExtTracksID>;
} // namespace o2::aod
