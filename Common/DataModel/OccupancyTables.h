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

/// \file OccupancyTables.h
/// \brief  Occupancy Table Header : TPC PID - Calibration
///
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (marian.ivanov@cern.ch)

#include <vector>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#ifndef COMMON_DATAMODEL_OCCUPANCYTABLES_H_
#define COMMON_DATAMODEL_OCCUPANCYTABLES_H_

namespace o2::aod
{
namespace occp
{
DECLARE_SOA_COLUMN(TfId, tfId, int64_t);
DECLARE_SOA_COLUMN(BcsInTFList, bcsInTFList, std::vector<int64_t>);

DECLARE_SOA_COLUMN(OccPrimUnfm80, occPrimUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccFV0AUnfm80, occFV0AUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccFV0CUnfm80, occFV0CUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccFT0AUnfm80, occFT0AUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccFT0CUnfm80, occFT0CUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccFDDAUnfm80, occFDDAUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccFDDCUnfm80, occFDDCUnfm80, std::vector<float>);

DECLARE_SOA_COLUMN(OccNTrackITSUnfm80, occNTrackITSUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackTPCUnfm80, occNTrackTPCUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackTRDUnfm80, occNTrackTRDUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackTOFUnfm80, occNTrackTOFUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackSizeUnfm80, occNTrackSizeUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackTPCAUnfm80, occNTrackTPCAUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackTPCCUnfm80, occNTrackTPCCUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackITSTPCUnfm80, occNTrackITSTPCUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackITSTPCAUnfm80, occNTrackITSTPCAUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccNTrackITSTPCCUnfm80, occNTrackITSTPCCUnfm80, std::vector<float>);

DECLARE_SOA_COLUMN(OccMultNTracksHasITSUnfm80, occMultNTracksHasITSUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccMultNTracksHasTPCUnfm80, occMultNTracksHasTPCUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccMultNTracksHasTOFUnfm80, occMultNTracksHasTOFUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccMultNTracksHasTRDUnfm80, occMultNTracksHasTRDUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccMultNTracksITSOnlyUnfm80, occMultNTracksITSOnlyUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccMultNTracksTPCOnlyUnfm80, occMultNTracksTPCOnlyUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccMultNTracksITSTPCUnfm80, occMultNTracksITSTPCUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccMultAllTracksTPCOnlyUnfm80, occMultAllTracksTPCOnlyUnfm80, std::vector<float>);

DECLARE_SOA_COLUMN(OccRobustT0V0PrimUnfm80, occRobustT0V0PrimUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccRobustFDDT0V0PrimUnfm80, occRobustFDDT0V0PrimUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccRobustNtrackDetUnfm80, occRobustNtrackDetUnfm80, std::vector<float>);
DECLARE_SOA_COLUMN(OccRobustMultExtraTableUnfm80, occRobustMultExtraTableUnfm80, std::vector<float>);

DECLARE_SOA_COLUMN(MeanOccPrimUnfm80, meanOccPrimUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFV0AUnfm80, meanOccFV0AUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFV0CUnfm80, meanOccFV0CUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFT0AUnfm80, meanOccFT0AUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFT0CUnfm80, meanOccFT0CUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFDDAUnfm80, meanOccFDDAUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFDDCUnfm80, meanOccFDDCUnfm80, float);

DECLARE_SOA_COLUMN(MeanOccNTrackITSUnfm80, meanOccNTrackITSUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTPCUnfm80, meanOccNTrackTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTRDUnfm80, meanOccNTrackTRDUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTOFUnfm80, meanOccNTrackTOFUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackSizeUnfm80, meanOccNTrackSizeUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTPCAUnfm80, meanOccNTrackTPCAUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTPCCUnfm80, meanOccNTrackTPCCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackITSTPCUnfm80, meanOccNTrackITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackITSTPCAUnfm80, meanOccNTrackITSTPCAUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackITSTPCCUnfm80, meanOccNTrackITSTPCCUnfm80, float);

DECLARE_SOA_COLUMN(MeanOccMultNTracksHasITSUnfm80, meanOccMultNTracksHasITSUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksHasTPCUnfm80, meanOccMultNTracksHasTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksHasTOFUnfm80, meanOccMultNTracksHasTOFUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksHasTRDUnfm80, meanOccMultNTracksHasTRDUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksITSOnlyUnfm80, meanOccMultNTracksITSOnlyUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksTPCOnlyUnfm80, meanOccMultNTracksTPCOnlyUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksITSTPCUnfm80, meanOccMultNTracksITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultAllTracksTPCOnlyUnfm80, meanOccMultAllTracksTPCOnlyUnfm80, float);

DECLARE_SOA_COLUMN(MeanOccRobustT0V0PrimUnfm80, meanOccRobustT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccRobustFDDT0V0PrimUnfm80, meanOccRobustFDDT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccRobustNtrackDetUnfm80, meanOccRobustNtrackDetUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccRobustMultExtraTableUnfm80, meanOccRobustMultExtraTableUnfm80, float);
} // namespace occp

DECLARE_SOA_TABLE(OccsBCsList, "AOD", "OCCSBCSLIST", o2::soa::Index<>, o2::aod::occp::TfId, o2::aod::occp::BcsInTFList);
DECLARE_SOA_TABLE(OccsDet, "AOD", "OCCSDET", o2::soa::Index<>,
                  o2::aod::occp::OccPrimUnfm80,
                  o2::aod::occp::OccFV0AUnfm80,
                  o2::aod::occp::OccFV0CUnfm80,
                  o2::aod::occp::OccFT0AUnfm80,
                  o2::aod::occp::OccFT0CUnfm80,
                  o2::aod::occp::OccFDDAUnfm80,
                  o2::aod::occp::OccFDDCUnfm80);

DECLARE_SOA_TABLE(OccsTrackMult, "AOD", "OCCSTRACKMULT", o2::soa::Index<>,
                  o2::aod::occp::OccNTrackITSUnfm80,
                  o2::aod::occp::OccNTrackTPCUnfm80,
                  o2::aod::occp::OccNTrackTRDUnfm80,
                  o2::aod::occp::OccNTrackTOFUnfm80,
                  o2::aod::occp::OccNTrackSizeUnfm80,
                  o2::aod::occp::OccNTrackTPCAUnfm80,
                  o2::aod::occp::OccNTrackTPCCUnfm80,
                  o2::aod::occp::OccNTrackITSTPCUnfm80,
                  o2::aod::occp::OccNTrackITSTPCAUnfm80,
                  o2::aod::occp::OccNTrackITSTPCCUnfm80);

DECLARE_SOA_TABLE(OccsMultExtra, "AOD", "OCCSMULTEXTRA", o2::soa::Index<>,
                  o2::aod::occp::OccMultNTracksHasITSUnfm80,
                  o2::aod::occp::OccMultNTracksHasTPCUnfm80,
                  o2::aod::occp::OccMultNTracksHasTOFUnfm80,
                  o2::aod::occp::OccMultNTracksHasTRDUnfm80,
                  o2::aod::occp::OccMultNTracksITSOnlyUnfm80,
                  o2::aod::occp::OccMultNTracksTPCOnlyUnfm80,
                  o2::aod::occp::OccMultNTracksITSTPCUnfm80,
                  o2::aod::occp::OccMultAllTracksTPCOnlyUnfm80);

DECLARE_SOA_TABLE(OccsRobust, "AOD", "OCCSROBUST", o2::soa::Index<>,
                  o2::aod::occp::OccRobustT0V0PrimUnfm80,
                  o2::aod::occp::OccRobustFDDT0V0PrimUnfm80,
                  o2::aod::occp::OccRobustNtrackDetUnfm80,
                  o2::aod::occp::OccRobustMultExtraTableUnfm80);

DECLARE_SOA_TABLE(OccsMeanDet, "AOD", "OCCSMEANDET", o2::soa::Index<>,
                  o2::aod::occp::MeanOccPrimUnfm80,
                  o2::aod::occp::MeanOccFV0AUnfm80,
                  o2::aod::occp::MeanOccFV0CUnfm80,
                  o2::aod::occp::MeanOccFT0AUnfm80,
                  o2::aod::occp::MeanOccFT0CUnfm80,
                  o2::aod::occp::MeanOccFDDAUnfm80,
                  o2::aod::occp::MeanOccFDDCUnfm80);

DECLARE_SOA_TABLE(OccsMeanTrkMult, "AOD", "OCCSMEANTRKMULT", o2::soa::Index<>,
                  o2::aod::occp::MeanOccNTrackITSUnfm80,
                  o2::aod::occp::MeanOccNTrackTPCUnfm80,
                  o2::aod::occp::MeanOccNTrackTRDUnfm80,
                  o2::aod::occp::MeanOccNTrackTOFUnfm80,
                  o2::aod::occp::MeanOccNTrackSizeUnfm80,
                  o2::aod::occp::MeanOccNTrackTPCAUnfm80,
                  o2::aod::occp::MeanOccNTrackTPCCUnfm80,
                  o2::aod::occp::MeanOccNTrackITSTPCUnfm80,
                  o2::aod::occp::MeanOccNTrackITSTPCAUnfm80,
                  o2::aod::occp::MeanOccNTrackITSTPCCUnfm80);

DECLARE_SOA_TABLE(OccsMnMultExtra, "AOD", "OCCSMNMULTEXTRA", o2::soa::Index<>,
                  o2::aod::occp::MeanOccMultNTracksHasITSUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTPCUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTOFUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTRDUnfm80,
                  o2::aod::occp::MeanOccMultNTracksITSOnlyUnfm80,
                  o2::aod::occp::MeanOccMultNTracksTPCOnlyUnfm80,
                  o2::aod::occp::MeanOccMultNTracksITSTPCUnfm80,
                  o2::aod::occp::MeanOccMultAllTracksTPCOnlyUnfm80);

DECLARE_SOA_TABLE(OccsMeanRobust, "AOD", "OCCSMEANROBUST", o2::soa::Index<>,
                  o2::aod::occp::MeanOccRobustT0V0PrimUnfm80,
                  o2::aod::occp::MeanOccRobustFDDT0V0PrimUnfm80,
                  o2::aod::occp::MeanOccRobustNtrackDetUnfm80,
                  o2::aod::occp::MeanOccRobustMultExtraTableUnfm80);

using Occs = aod::OccsBCsList;
using Occ = Occs::iterator;

namespace occidx
{
DECLARE_SOA_INDEX_COLUMN(BC, bc); // Iterator is passed here in index column
DECLARE_SOA_INDEX_COLUMN(Occ, occ);
DECLARE_SOA_COLUMN(TfId, tfId, int);
DECLARE_SOA_COLUMN(BcInTF, bcInTF, int);
} // namespace occidx

// DECLARE_SOA_TABLE(OccIndexTable, "AOD", "OCCINDEXTABLE", o2::soa::Index<>,
DECLARE_SOA_INDEX_TABLE_USER(OccIndexTable, Occs, "OCCINDEXTABLE",
                             o2::aod::occidx::BCId,
                             o2::aod::occidx::OccId);

DECLARE_SOA_TABLE(BCTFinfoTable, "AOD", "BCTFINFOTABLE", o2::soa::Index<>,
                  o2::aod::occidx::TfId,
                  o2::aod::occidx::BcInTF);

namespace trackmeanocc
{
DECLARE_SOA_INDEX_COLUMN(Track, track);

DECLARE_SOA_COLUMN(MeanOccPrimUnfm80, meanOccPrimUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFV0AUnfm80, meanOccFV0AUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFV0CUnfm80, meanOccFV0CUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFT0AUnfm80, meanOccFT0AUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFT0CUnfm80, meanOccFT0CUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFDDAUnfm80, meanOccFDDAUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccFDDCUnfm80, meanOccFDDCUnfm80, float);

DECLARE_SOA_COLUMN(MeanOccNTrackITSUnfm80, meanOccNTrackITSUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTPCUnfm80, meanOccNTrackTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTRDUnfm80, meanOccNTrackTRDUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTOFUnfm80, meanOccNTrackTOFUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackSizeUnfm80, meanOccNTrackSizeUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTPCAUnfm80, meanOccNTrackTPCAUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackTPCCUnfm80, meanOccNTrackTPCCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackITSTPCUnfm80, meanOccNTrackITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackITSTPCAUnfm80, meanOccNTrackITSTPCAUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccNTrackITSTPCCUnfm80, meanOccNTrackITSTPCCUnfm80, float);

DECLARE_SOA_COLUMN(MeanOccMultNTracksHasITSUnfm80, meanOccMultNTracksHasITSUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksHasTPCUnfm80, meanOccMultNTracksHasTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksHasTOFUnfm80, meanOccMultNTracksHasTOFUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksHasTRDUnfm80, meanOccMultNTracksHasTRDUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksITSOnlyUnfm80, meanOccMultNTracksITSOnlyUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksTPCOnlyUnfm80, meanOccMultNTracksTPCOnlyUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultNTracksITSTPCUnfm80, meanOccMultNTracksITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccMultAllTracksTPCOnlyUnfm80, meanOccMultAllTracksTPCOnlyUnfm80, float);

DECLARE_SOA_COLUMN(MeanOccRobustT0V0PrimUnfm80, meanOccRobustT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccRobustFDDT0V0PrimUnfm80, meanOccRobustFDDT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccRobustNtrackDetUnfm80, meanOccRobustNtrackDetUnfm80, float);
DECLARE_SOA_COLUMN(MeanOccRobustMultExtraTableUnfm80, meanOccRobustMultExtraTableUnfm80, float);

DECLARE_SOA_COLUMN(WeightMeanOccPrimUnfm80, weightMeanOccPrimUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccFV0AUnfm80, weightMeanOccFV0AUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccFV0CUnfm80, weightMeanOccFV0CUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccFT0AUnfm80, weightMeanOccFT0AUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccFT0CUnfm80, weightMeanOccFT0CUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccFDDAUnfm80, weightMeanOccFDDAUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccFDDCUnfm80, weightMeanOccFDDCUnfm80, float);

DECLARE_SOA_COLUMN(WeightMeanOccNTrackITSUnfm80, weightMeanOccNTrackITSUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackTPCUnfm80, weightMeanOccNTrackTPCUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackTRDUnfm80, weightMeanOccNTrackTRDUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackTOFUnfm80, weightMeanOccNTrackTOFUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackSizeUnfm80, weightMeanOccNTrackSizeUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackTPCAUnfm80, weightMeanOccNTrackTPCAUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackTPCCUnfm80, weightMeanOccNTrackTPCCUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackITSTPCUnfm80, weightMeanOccNTrackITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackITSTPCAUnfm80, weightMeanOccNTrackITSTPCAUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccNTrackITSTPCCUnfm80, weightMeanOccNTrackITSTPCCUnfm80, float);

DECLARE_SOA_COLUMN(WeightMeanOccMultNTracksHasITSUnfm80, weightMeanOccMultNTracksHasITSUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccMultNTracksHasTPCUnfm80, weightMeanOccMultNTracksHasTPCUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccMultNTracksHasTOFUnfm80, weightMeanOccMultNTracksHasTOFUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccMultNTracksHasTRDUnfm80, weightMeanOccMultNTracksHasTRDUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccMultNTracksITSOnlyUnfm80, weightMeanOccMultNTracksITSOnlyUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccMultNTracksTPCOnlyUnfm80, weightMeanOccMultNTracksTPCOnlyUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccMultNTracksITSTPCUnfm80, weightMeanOccMultNTracksITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccMultAllTracksTPCOnlyUnfm80, weightMeanOccMultAllTracksTPCOnlyUnfm80, float);

DECLARE_SOA_COLUMN(WeightMeanOccRobustT0V0PrimUnfm80, weightMeanOccRobustT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccRobustFDDT0V0PrimUnfm80, weightMeanOccRobustFDDT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccRobustNtrackDetUnfm80, weightMeanOccRobustNtrackDetUnfm80, float);
DECLARE_SOA_COLUMN(WeightMeanOccRobustMultExtraTableUnfm80, weightMeanOccRobustMultExtraTableUnfm80, float);

} // namespace trackmeanocc

// Tracks
// using Tracks = aod::Tracks;
// DECLARE_SOA_INDEX_TABLE_USER(TrackMeanOccs0, Tracks, "TRACKMEANOCCS0", o2::aod::trackmeanocc::TrackId);
DECLARE_SOA_TABLE(TrackMeanOccs0, "AOD", "TRACKMEANOCCS0", o2::aod::trackmeanocc::TrackId);

DECLARE_SOA_TABLE(TrackMeanOccs1, "AOD", "TRACKMEANOCCS1", o2::soa::Index<>, // TrackMeanOccDet
                  o2::aod::trackmeanocc::MeanOccPrimUnfm80,
                  o2::aod::trackmeanocc::MeanOccFV0AUnfm80,
                  o2::aod::trackmeanocc::MeanOccFV0CUnfm80,
                  o2::aod::trackmeanocc::MeanOccFT0AUnfm80,
                  o2::aod::trackmeanocc::MeanOccFT0CUnfm80,
                  o2::aod::trackmeanocc::MeanOccFDDAUnfm80,
                  o2::aod::trackmeanocc::MeanOccFDDCUnfm80);

DECLARE_SOA_TABLE(TrackMeanOccs2, "AOD", "TRACKMEANOCCS2", o2::soa::Index<>, // TrackMeanOccTrackMult
                  o2::aod::trackmeanocc::MeanOccNTrackITSUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::MeanOccNTrackITSTPCCUnfm80);

DECLARE_SOA_TABLE(TrackMeanOccs3, "AOD", "TRACKMEANOCCS3", o2::soa::Index<>, // TrackMeanOccMultExtra
                  o2::aod::trackmeanocc::MeanOccMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::MeanOccMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::MeanOccMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::MeanOccMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::MeanOccMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::MeanOccMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::MeanOccMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::MeanOccMultAllTracksTPCOnlyUnfm80);

DECLARE_SOA_TABLE(TrackMeanOccs4, "AOD", "TRACKMEANOCCS4", o2::soa::Index<>, // TrackMeanOccRobus
                  o2::aod::trackmeanocc::MeanOccRobustT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::MeanOccRobustFDDT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::MeanOccRobustNtrackDetUnfm80,
                  o2::aod::trackmeanocc::MeanOccRobustMultExtraTableUnfm80);

DECLARE_SOA_TABLE(TrackMeanOccs5, "AOD", "TRACKMEANOCCS5", o2::soa::Index<>, // TrackWieghtMeanOccDet
                  o2::aod::trackmeanocc::WeightMeanOccPrimUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccFV0AUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccFV0CUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccFT0AUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccFT0CUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccFDDAUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccFDDCUnfm80);

DECLARE_SOA_TABLE(TrackMeanOccs6, "AOD", "TRACKMEANOCCS6", o2::soa::Index<>, // TrackWieghtMeanOccMult
                  o2::aod::trackmeanocc::WeightMeanOccNTrackITSUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccNTrackITSTPCCUnfm80)

DECLARE_SOA_TABLE(TrackMeanOccs7, "AOD", "TRACKMEANOCCS7", o2::soa::Index<>, // TrackWeightMeanOccMultExtra
                  o2::aod::trackmeanocc::WeightMeanOccMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccMultAllTracksTPCOnlyUnfm80);

DECLARE_SOA_TABLE(TrackMeanOccs8, "AOD", "TRACKMEANOCCS8", o2::soa::Index<>, // TrackWieghtMeanOccRboust
                  o2::aod::trackmeanocc::WeightMeanOccRobustT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccRobustFDDT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccRobustNtrackDetUnfm80,
                  o2::aod::trackmeanocc::WeightMeanOccRobustMultExtraTableUnfm80);
} // namespace o2::aod
#endif // COMMON_DATAMODEL_OCCUPANCYTABLES_H_
