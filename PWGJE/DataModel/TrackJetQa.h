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
/// \file   JetTrackQa.h
/// \author Johanna Lömker <johanna.lomker@cern.ch>
/// \since  2023-10-02
/// \brief  Header for the trackJetQa task for the analysis of the tracks for jets.
///

#ifndef PWGJE_DATAMODEL_TRACKJETQA_H_
#define PWGJE_DATAMODEL_TRACKJETQA_H_

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/Jet.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Derived data model for track optimization (and cut variation)
namespace o2::aod
{
namespace jetcollisions
{
// Collision info
DECLARE_SOA_COLUMN(GlobalIdx, globalIdx, int64_t);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(MultNTracksPV, multNTracksPV, int);
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float);
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);
} // namespace jetcollisions

namespace jettrack
{
// Track info
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);          //! Id of collision
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool); //! IsPVContributor
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);                   //! Has or not the TRD match
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);                   //! Has or not the ITS match
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);                   //! Has or not the TPC match         // Addtional selections for trackJetqa //
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);                                   // if a track passed the isGlobalTrack requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);                         // if a track passed the isGlobalTrackWoDCA requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoPtEta, isGlobalTrackWoPtEta, bool);                     // if a track passed the isGlobalTrackWoPtEta requirement
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);
DECLARE_SOA_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, float);
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, float);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int16_t);
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(TPCFoundOverFindableCls, tpcFoundOverFindableCls, float);

} // namespace jettrack

DECLARE_SOA_TABLE(JeColls, "AOD", "JECOLLS",
                  o2::soa::Index<>,
                  collision::CollisionTime,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  jetcollisions::Sel8,
                  jetcollisions::MultNTracksPV,
                  jetcollisions::MultFT0A,
                  jetcollisions::MultFT0C,
                  jetcollisions::CentFT0A,
                  jetcollisions::CentFT0C,
                  jetcollisions::RunNumber);

// using JeColl = JeColls::iterator;

DECLARE_SOA_TABLE(JeTracks, "AOD", "JETRACKS",
                  o2::soa::Index<>,
                  jettrack::CollisionId,
                  track::TrackTime,
                  track::Signed1Pt, track::Eta, track::Phi, track::Pt,
                  track::Sigma1Pt,
                  track::Alpha,
                  track::X, track::Y, track::Z,
                  track::Snp,
                  track::Tgl,
                  jettrack::IsPVContributor,
                  jettrack::HasTRD,
                  jettrack::HasITS,
                  jettrack::HasTPC,
                  jettrack::IsGlobalTrack,
                  jettrack::IsGlobalTrackWoDCA,
                  jettrack::IsGlobalTrackWoPtEta,
                  track::Flags,
                  track::TrackType,
                  track::Length,
                  track::TPCChi2NCl, track::ITSChi2NCl, track::TOFChi2,
                  track::TPCNClsShared,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  track::ITSClusterMap,
                  jettrack::ITSNCls,
                  jettrack::TPCFractionSharedCls,
                  jettrack::TPCNClsFound,
                  jettrack::TPCNClsCrossedRows,
                  jettrack::TPCCrossedRowsOverFindableCls,
                  jettrack::TPCFoundOverFindableCls,
                  track::DcaXY,
                  track::DcaZ);
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_TRACKJETQA_H_
