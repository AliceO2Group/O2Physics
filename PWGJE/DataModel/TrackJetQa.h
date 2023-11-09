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
namespace jetspectra
{

// Collision info
DECLARE_SOA_INDEX_COLUMN(BC, bc); //! Most probably BC to where this collision has occurred
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(MultNTracksPV, multNTracksPV, int);
DECLARE_SOA_COLUMN(MultTracklets, multTracklets, int);
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
// Track info
DECLARE_SOA_INDEX_COLUMN(Collision, collision);             //! Index to the collision
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool); //! IsPVContributor
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);                   //! Has or not the TRD match
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);                                   // if a track passed the isGlobalTrack requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);                         // if a track passed the isGlobalTrackWoDCA requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoPtEta, isGlobalTrackWoPtEta, bool);                     // all tracks in the derived table have to pass this requirement !
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);
DECLARE_SOA_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, float);
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, float);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int16_t);
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(TPCFoundOverFindableCls, tpcFoundOverFindableCls, float);

} // namespace jetspectra

DECLARE_SOA_TABLE(JeColls, "AOD", "JECOLLS",
                  o2::soa::Index<>,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  jetspectra::Sel8,
                  jetspectra::MultNTracksPV,
                  jetspectra::MultTracklets,
                  jetspectra::MultFT0M,
                  jetspectra::CentFT0M,
                  jetspectra::RunNumber);
using JeColl = JeColls::iterator;

DECLARE_SOA_TABLE(JeTracks, "AOD", "JETRACKS",
                  o2::soa::Index<>,
                  jetspectra::CollisionId,
                  track::Signed1Pt, track::Eta, track::Phi, track::Pt,
                  track::Sigma1Pt,
                  track::Alpha,
                  track::X, track::Y, track::Z,
                  track::Snp,
                  track::Tgl,
                  jetspectra::IsPVContributor,
                  jetspectra::HasTRD,
                  jetspectra::IsGlobalTrack,
                  jetspectra::IsGlobalTrackWoDCA,
                  jetspectra::IsGlobalTrackWoPtEta,
                  track::Flags,
                  track::Length,
                  track::TPCChi2NCl, track::ITSChi2NCl, track::TOFChi2,
                  track::TPCNClsShared,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  track::ITSClusterMap,
                  jetspectra::ITSNCls,
                  jetspectra::TPCFractionSharedCls,
                  jetspectra::TPCNClsFound,
                  jetspectra::TPCNClsCrossedRows,
                  jetspectra::TPCCrossedRowsOverFindableCls,
                  jetspectra::TPCFoundOverFindableCls,
                  track::DcaXY,
                  track::DcaZ);
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_TRACKJETQA_H_
