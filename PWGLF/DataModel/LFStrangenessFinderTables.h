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
//
// This file defines the auxiliary tables to be used in the v0 and cascade
// finders. These are cross-check tasks that are not meant to do final analyses
// as finding will be extremely slow and complex at the AO2D level.

#ifndef PWGLF_DATAMODEL_LFSTRANGENESSFINDERTABLES_H_
#define PWGLF_DATAMODEL_LFSTRANGENESSFINDERTABLES_H_

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>

// V0 auxiliary tables
namespace o2::aod
{
namespace vFinderTrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Track");
DECLARE_SOA_COLUMN(IsPositive, isPositive, bool);     // is positively charged
DECLARE_SOA_COLUMN(CompatiblePi, compatiblePi, bool); // compatible with pion hypothesis
DECLARE_SOA_COLUMN(CompatibleKa, compatibleKa, bool); // compatible with kaon hypothesis
DECLARE_SOA_COLUMN(CompatiblePr, compatiblePr, bool); // compatible with proton hypothesis
} // namespace vFinderTrack
DECLARE_SOA_TABLE(VFinderTracks, "AOD", "VFINDERTRACK", o2::soa::Index<>, vFinderTrack::TrackId, vFinderTrack::IsPositive, vFinderTrack::CompatiblePi, vFinderTrack::CompatibleKa, vFinderTrack::CompatiblePr);

// Cascade auxiliary tables
// Store good bachelor track candidates
// namespace cFinderTrack
// {
// DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Track");
// DECLARE_SOA_COLUMN(IsPositive, isPositive, bool); // is positively charged
// DECLARE_SOA_COLUMN(CompatiblePi, compatiblePi, bool); // compatible with pion hypothesis
// DECLARE_SOA_COLUMN(CompatibleKa, compatibleKa, bool); // compatible with pion hypothesis
// DECLARE_SOA_COLUMN(CompatiblePr, compatiblePr, bool); // compatible with pion hypothesis
// } // namespace cFinderTrack
// DECLARE_SOA_TABLE(CFinderTracks, "AOD", "CFINDERTRACK", o2::soa::Index<>, cFinderTrack::TrackId, cFinderTrack::IsPositive, cFinderTrack::CompatiblePi, cFinderTrack::CompatibleKa, cFinderTrack::CompatiblePr);

// // Store good V0 candidates
// namespace cFinderV0
// {
// DECLARE_SOA_INDEX_COLUMN_FULL(Lambda, lambda, int, V0Datas, "_Lambda");
// DECLARE_SOA_INDEX_COLUMN(CompatibleLamba, compatibleLambda); // compatible with Lambda
// DECLARE_SOA_INDEX_COLUMN(CompatibleAntiLamba, compatibleAntiLambda); // compatible with AntiLambda
// } // namespace cFinderV0
// DECLARE_SOA_TABLE(CFinderV0, "AOD", "VFINDERV0", o2::soa::Index<>, cFinderV0::LambdaId, cFinderV0::CompatibleLamba, cFinderV0::CompatibleAntiLamba);

namespace cascgoodpostracks
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodPosTrack, goodPosTrack, int, Tracks, "_GoodPos");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
} // namespace cascgoodpostracks
DECLARE_SOA_TABLE(CascGoodPosTracks, "AOD", "CASCGOODPTRACKS", o2::soa::Index<>, cascgoodpostracks::GoodPosTrackId, cascgoodpostracks::CollisionId, cascgoodpostracks::DCAXY);
namespace cascgoodnegtracks
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodNegTrack, goodNegTrack, int, Tracks, "_GoodNeg");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
} // namespace cascgoodnegtracks
DECLARE_SOA_TABLE(CascGoodNegTracks, "AOD", "CASCGOODNTRACKS", o2::soa::Index<>, cascgoodnegtracks::GoodNegTrackId, cascgoodnegtracks::CollisionId, cascgoodnegtracks::DCAXY);
namespace cascgoodlambdas
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodLambda, goodLambda, int, V0Datas, "_GoodLambda");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace cascgoodlambdas
DECLARE_SOA_TABLE(CascGoodLambdas, "AOD", "CASCGOODLAM", o2::soa::Index<>, cascgoodlambdas::GoodLambdaId, cascgoodlambdas::CollisionId);
namespace cascgoodantilambdas
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodAntiLambda, goodAntiLambda, int, V0Datas, "_GoodAntiLambda");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace cascgoodantilambdas
DECLARE_SOA_TABLE(CascGoodAntiLambdas, "AOD", "CASCGOODALAM", o2::soa::Index<>, cascgoodantilambdas::GoodAntiLambdaId, cascgoodantilambdas::CollisionId);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSTRANGENESSFINDERTABLES_H_
