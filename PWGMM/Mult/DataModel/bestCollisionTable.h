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
///
/// \file   bestCollisionTable.h
/// \brief This code produces tables including central and MFT tracks based on smallest DCAxy/DCAz approach
/// \author Anton Alkin <anton.alkin@cern.ch>
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>
/// \author Gyula Bencedi <gyula.bencedi@cern.ch>
/// \author Tulika Tripathy <tulika.tripathy@cern.ch>

#ifndef PWGMM_MULT_DATAMODEL_BESTCOLLISIONTABLE_H_
#define PWGMM_MULT_DATAMODEL_BESTCOLLISIONTABLE_H_

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
DECLARE_SOA_INDEX_COLUMN_FULL(BestCollision, bestCollision, int, Collisions, "");
DECLARE_SOA_COLUMN(AmbDegree, ambDegree, int); // degree of ambiguity of the track
DECLARE_SOA_COLUMN(BestDCAXY, bestDCAXY, float);
DECLARE_SOA_COLUMN(BestDCAX, bestDCAX, float);
DECLARE_SOA_COLUMN(BestDCAY, bestDCAY, float);
DECLARE_SOA_COLUMN(BestDCAZ, bestDCAZ, float);
DECLARE_SOA_COLUMN(PtStatic, pts, float);
DECLARE_SOA_COLUMN(PStatic, ps, float);
DECLARE_SOA_COLUMN(EtaStatic, etas, float);
DECLARE_SOA_COLUMN(PhiStatic, phis, float);
} // namespace fwdtrack

namespace pwgmm::indices
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(MFTTrack, mfttrack);
} // namespace pwgmm::indices

DECLARE_SOA_TABLE(BestCollisionsFwd, "AOD", "BESTCOLLFWD", o2::soa::Index<>, pwgmm::indices::MFTTrackId, aod::fwdtrack::AmbDegree,
                  aod::fwdtrack::BestCollisionId, aod::fwdtrack::BestDCAXY,
                  fwdtrack::BestDCAX, fwdtrack::BestDCAY); // beware: depending on which process produced this table,
// it can be joined with either MFTAmbiguousTracks OR MFTTracks

DECLARE_SOA_TABLE(BestCollFwdExtra, "AOD", "BESTCOLLFWDE",
                  fwdtrack::X, fwdtrack::Y,
                  fwdtrack::Z, fwdtrack::Tgl, fwdtrack::Signed1Pt,
                  fwdtrack::PtStatic, fwdtrack::PStatic, fwdtrack::EtaStatic,
                  fwdtrack::PhiStatic); // Snp does not exist

DECLARE_SOA_TABLE(BestCollisionsFwd3d, "AOD", "BESTCOLLFWD3D",
                  o2::soa::Index<>,
                  pwgmm::indices::MFTTrackId,
                  aod::fwdtrack::AmbDegree,
                  aod::fwdtrack::BestCollisionId,
                  aod::fwdtrack::BestDCAXY,
                  aod::fwdtrack::BestDCAZ);

DECLARE_SOA_TABLE(BestCollisionsFwd3dExtra, "AOD", "BESTCOLLFWD3DE",
                  fwdtrack::X, fwdtrack::Y,
                  fwdtrack::Z, fwdtrack::Tgl, fwdtrack::Signed1Pt,
                  fwdtrack::PtStatic, fwdtrack::PStatic, fwdtrack::EtaStatic,
                  fwdtrack::PhiStatic); // Snp does not exist

DECLARE_SOA_TABLE(ReassignedTracksCore, "AOD", "CRRETRACKS",
                  aod::track::BestCollisionId,
                  pwgmm::indices::TrackId,
                  aod::track::BestDCAXY,
                  aod::track::BestDCAZ);

DECLARE_SOA_TABLE(ReassignedTracksExtra, "AOD", "EXRETRACKS",
                  track::X, track::Alpha, track::Y,
                  track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::PtStatic, track::PStatic, track::EtaStatic,
                  track::PhiStatic);

} // namespace o2::aod

#endif // PWGMM_MULT_DATAMODEL_BESTCOLLISIONTABLE_H_
