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

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace track
{
DECLARE_SOA_INDEX_COLUMN_FULL(BestCollision, bestCollision, int32_t, Collisions, "");
DECLARE_SOA_COLUMN(BestDCAXY, bestDCAXY, float);
DECLARE_SOA_COLUMN(BestDCAZ, bestDCAZ, float);
DECLARE_SOA_COLUMN(PtStatic, pts, float);
DECLARE_SOA_COLUMN(PStatic, ps, float);
DECLARE_SOA_COLUMN(EtaStatic, etas, float);
DECLARE_SOA_COLUMN(PhiStatic, phis, float);
} // namespace track
namespace fwdtrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(BestCollision, bestCollision, int32_t, Collisions, "");
DECLARE_SOA_COLUMN(BestDCAXY, bestDCAXY, float);
DECLARE_SOA_COLUMN(PtStatic, pts, float);
DECLARE_SOA_COLUMN(PStatic, ps, float);
DECLARE_SOA_COLUMN(EtaStatic, etas, float);
DECLARE_SOA_COLUMN(PhiStatic, phis, float);
} // namespace fwdtrack
DECLARE_SOA_TABLE(BestCollisions, "AOD", "BESTCOLL",
                  aod::track::BestCollisionId, aod::track::BestDCAXY,
                  aod::track::BestDCAZ, track::X, track::Alpha, track::Y,
                  track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::PtStatic, track::PStatic, track::EtaStatic,
                  track::PhiStatic);

DECLARE_SOA_TABLE(BestCollisionsFwd, "AOD", "BESTCOLLFWD",
                  aod::fwdtrack::BestCollisionId, aod::fwdtrack::BestDCAXY,
                  fwdtrack::X, fwdtrack::Y,
                  fwdtrack::Z, fwdtrack::Tgl, fwdtrack::Signed1Pt,
                  fwdtrack::PtStatic, fwdtrack::PStatic, fwdtrack::EtaStatic,
                  fwdtrack::PhiStatic); // Snp does not exist
} // namespace o2::aod
