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

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>
#include <vector>

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
// 1
DECLARE_SOA_TABLE(OccsPrim, "AOD", "OCCSPRIM", o2::soa::Index<>,
                  o2::aod::occp::OccPrimUnfm80);
DECLARE_SOA_TABLE(OccsMeanPrim, "AOD", "OCCSMEANPRIM", o2::soa::Index<>,
                  o2::aod::occp::MeanOccPrimUnfm80);
// 2
DECLARE_SOA_TABLE(OccsT0V0, "AOD", "OCCST0V0", o2::soa::Index<>,
                  o2::aod::occp::OccFV0AUnfm80,
                  o2::aod::occp::OccFV0CUnfm80,
                  o2::aod::occp::OccFT0AUnfm80,
                  o2::aod::occp::OccFT0CUnfm80);
DECLARE_SOA_TABLE(OccsMeanT0V0, "AOD", "OCCSMEANT0V0", o2::soa::Index<>,
                  o2::aod::occp::MeanOccFV0AUnfm80,
                  o2::aod::occp::MeanOccFV0CUnfm80,
                  o2::aod::occp::MeanOccFT0AUnfm80,
                  o2::aod::occp::MeanOccFT0CUnfm80);
// 3
DECLARE_SOA_TABLE(OccsFDD, "AOD", "OCCSFDD", o2::soa::Index<>,
                  o2::aod::occp::OccFDDAUnfm80,
                  o2::aod::occp::OccFDDCUnfm80);
DECLARE_SOA_TABLE(OccsMeanFDD, "AOD", "OCCSMEANFDD", o2::soa::Index<>,
                  o2::aod::occp::MeanOccFDDAUnfm80,
                  o2::aod::occp::MeanOccFDDCUnfm80);
// 4
DECLARE_SOA_TABLE(OccsNTrackDet, "AOD", "OCCSNTRACKDET", o2::soa::Index<>,
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
DECLARE_SOA_TABLE(OccsMeanNTrkDet, "AOD", "OCCSMEANNTRKDET", o2::soa::Index<>,
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
// 5
DECLARE_SOA_TABLE(OccsMultExtra, "AOD", "OCCSMULTEXTRA", o2::soa::Index<>,
                  o2::aod::occp::OccMultNTracksHasITSUnfm80,
                  o2::aod::occp::OccMultNTracksHasTPCUnfm80,
                  o2::aod::occp::OccMultNTracksHasTOFUnfm80,
                  o2::aod::occp::OccMultNTracksHasTRDUnfm80,
                  o2::aod::occp::OccMultNTracksITSOnlyUnfm80,
                  o2::aod::occp::OccMultNTracksTPCOnlyUnfm80,
                  o2::aod::occp::OccMultNTracksITSTPCUnfm80,
                  o2::aod::occp::OccMultAllTracksTPCOnlyUnfm80);
DECLARE_SOA_TABLE(OccsMnMultExtra, "AOD", "OCCSMNMULTEXTRA", o2::soa::Index<>,
                  o2::aod::occp::MeanOccMultNTracksHasITSUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTPCUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTOFUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTRDUnfm80,
                  o2::aod::occp::MeanOccMultNTracksITSOnlyUnfm80,
                  o2::aod::occp::MeanOccMultNTracksTPCOnlyUnfm80,
                  o2::aod::occp::MeanOccMultNTracksITSTPCUnfm80,
                  o2::aod::occp::MeanOccMultAllTracksTPCOnlyUnfm80);

// Robust Occupancies
// 1
DECLARE_SOA_TABLE(ORT0V0Prim, "AOD", "ORT0V0Prim", o2::soa::Index<>,
                  o2::aod::occp::OccRobustT0V0PrimUnfm80);
DECLARE_SOA_TABLE(OMRT0V0Prim, "AOD", "OMRT0V0PRIM", o2::soa::Index<>,
                  o2::aod::occp::MeanOccRobustT0V0PrimUnfm80);
// 2
DECLARE_SOA_TABLE(ORFDDT0V0Prim, "AOD", "ORFDDT0V0PRIM", o2::soa::Index<>,
                  o2::aod::occp::OccRobustFDDT0V0PrimUnfm80);
DECLARE_SOA_TABLE(OMRFDDT0V0Prim, "AOD", "OMRFDDT0V0PRIM", o2::soa::Index<>,
                  o2::aod::occp::MeanOccRobustFDDT0V0PrimUnfm80);
// 3
DECLARE_SOA_TABLE(ORNtrackDet, "AOD", "ORNTRACKDET", o2::soa::Index<>,
                  o2::aod::occp::OccRobustNtrackDetUnfm80);
DECLARE_SOA_TABLE(OMRNtrackDet, "AOD", "OMRNTRACKDET", o2::soa::Index<>,
                  o2::aod::occp::MeanOccRobustNtrackDetUnfm80);
// 4
DECLARE_SOA_TABLE(ORMultExtra, "AOD", "ORMULTEXTRA", o2::soa::Index<>,
                  o2::aod::occp::OccRobustMultExtraTableUnfm80);
DECLARE_SOA_TABLE(OMRMultExtra, "AOD", "OMRMULTEXTRA", o2::soa::Index<>,
                  o2::aod::occp::MeanOccRobustMultExtraTableUnfm80);

using Occs = aod::OccsBCsList;
using Occ = Occs::iterator;

namespace occidx
{
DECLARE_SOA_INDEX_COLUMN(BC, bc); // o2-linter: disable=name/o2-column (BC is an acronym already defined in data model) // Iterator is passed here in index column
DECLARE_SOA_INDEX_COLUMN(Occ, occ);
DECLARE_SOA_COLUMN(TfId, tfId, int);
DECLARE_SOA_COLUMN(BcInTF, bcInTF, int);
} // namespace occidx

DECLARE_SOA_TABLE(OccIndexTable, "AOD", "OCCINDEXTABLE", o2::soa::Index<>,
                  o2::aod::occidx::BCId,
                  o2::aod::occidx::OccId);

DECLARE_SOA_TABLE(BCTFinfoTable, "AOD", "BCTFINFOTABLE", o2::soa::Index<>,
                  o2::aod::occidx::TfId,
                  o2::aod::occidx::BcInTF);

namespace trackmeanocc
{
DECLARE_SOA_INDEX_COLUMN(Track, track);

// DECLARE_SOA_INDEX_COLUMN(TracksQA, tracksQA);// this is not working

DECLARE_SOA_COLUMN(TmoPrimUnfm80, tmoPrimUnfm80, float);
DECLARE_SOA_COLUMN(TmoFV0AUnfm80, tmoFV0AUnfm80, float);
DECLARE_SOA_COLUMN(TmoFV0CUnfm80, tmoFV0CUnfm80, float);
DECLARE_SOA_COLUMN(TmoFT0AUnfm80, tmoFT0AUnfm80, float);
DECLARE_SOA_COLUMN(TmoFT0CUnfm80, tmoFT0CUnfm80, float);
DECLARE_SOA_COLUMN(TmoFDDAUnfm80, tmoFDDAUnfm80, float);
DECLARE_SOA_COLUMN(TmoFDDCUnfm80, tmoFDDCUnfm80, float);

DECLARE_SOA_COLUMN(TmoNTrackITSUnfm80, tmoNTrackITSUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackTPCUnfm80, tmoNTrackTPCUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackTRDUnfm80, tmoNTrackTRDUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackTOFUnfm80, tmoNTrackTOFUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackSizeUnfm80, tmoNTrackSizeUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackTPCAUnfm80, tmoNTrackTPCAUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackTPCCUnfm80, tmoNTrackTPCCUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackITSTPCUnfm80, tmoNTrackITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackITSTPCAUnfm80, tmoNTrackITSTPCAUnfm80, float);
DECLARE_SOA_COLUMN(TmoNTrackITSTPCCUnfm80, tmoNTrackITSTPCCUnfm80, float);

DECLARE_SOA_COLUMN(TmoMultNTracksHasITSUnfm80, tmoMultNTracksHasITSUnfm80, float);
DECLARE_SOA_COLUMN(TmoMultNTracksHasTPCUnfm80, tmoMultNTracksHasTPCUnfm80, float);
DECLARE_SOA_COLUMN(TmoMultNTracksHasTOFUnfm80, tmoMultNTracksHasTOFUnfm80, float);
DECLARE_SOA_COLUMN(TmoMultNTracksHasTRDUnfm80, tmoMultNTracksHasTRDUnfm80, float);
DECLARE_SOA_COLUMN(TmoMultNTracksITSOnlyUnfm80, tmoMultNTracksITSOnlyUnfm80, float);
DECLARE_SOA_COLUMN(TmoMultNTracksTPCOnlyUnfm80, tmoMultNTracksTPCOnlyUnfm80, float);
DECLARE_SOA_COLUMN(TmoMultNTracksITSTPCUnfm80, tmoMultNTracksITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(TmoMultAllTracksTPCOnlyUnfm80, tmoMultAllTracksTPCOnlyUnfm80, float);

DECLARE_SOA_COLUMN(TmoRobustT0V0PrimUnfm80, tmoRobustT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(TmoRobustFDDT0V0PrimUnfm80, tmoRobustFDDT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(TmoRobustNtrackDetUnfm80, tmoRobustNtrackDetUnfm80, float);
DECLARE_SOA_COLUMN(TmoRobustMultExtraTableUnfm80, tmoRobustMultExtraTableUnfm80, float);

DECLARE_SOA_COLUMN(TwmoPrimUnfm80, twmoPrimUnfm80, float);
DECLARE_SOA_COLUMN(TwmoFV0AUnfm80, twmoFV0AUnfm80, float);
DECLARE_SOA_COLUMN(TwmoFV0CUnfm80, twmoFV0CUnfm80, float);
DECLARE_SOA_COLUMN(TwmoFT0AUnfm80, twmoFT0AUnfm80, float);
DECLARE_SOA_COLUMN(TwmoFT0CUnfm80, twmoFT0CUnfm80, float);
DECLARE_SOA_COLUMN(TwmoFDDAUnfm80, twmoFDDAUnfm80, float);
DECLARE_SOA_COLUMN(TwmoFDDCUnfm80, twmoFDDCUnfm80, float);

DECLARE_SOA_COLUMN(TwmoNTrackITSUnfm80, twmoNTrackITSUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackTPCUnfm80, twmoNTrackTPCUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackTRDUnfm80, twmoNTrackTRDUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackTOFUnfm80, twmoNTrackTOFUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackSizeUnfm80, twmoNTrackSizeUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackTPCAUnfm80, twmoNTrackTPCAUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackTPCCUnfm80, twmoNTrackTPCCUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackITSTPCUnfm80, twmoNTrackITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackITSTPCAUnfm80, twmoNTrackITSTPCAUnfm80, float);
DECLARE_SOA_COLUMN(TwmoNTrackITSTPCCUnfm80, twmoNTrackITSTPCCUnfm80, float);

DECLARE_SOA_COLUMN(TwmoMultNTracksHasITSUnfm80, twmoMultNTracksHasITSUnfm80, float);
DECLARE_SOA_COLUMN(TwmoMultNTracksHasTPCUnfm80, twmoMultNTracksHasTPCUnfm80, float);
DECLARE_SOA_COLUMN(TwmoMultNTracksHasTOFUnfm80, twmoMultNTracksHasTOFUnfm80, float);
DECLARE_SOA_COLUMN(TwmoMultNTracksHasTRDUnfm80, twmoMultNTracksHasTRDUnfm80, float);
DECLARE_SOA_COLUMN(TwmoMultNTracksITSOnlyUnfm80, twmoMultNTracksITSOnlyUnfm80, float);
DECLARE_SOA_COLUMN(TwmoMultNTracksTPCOnlyUnfm80, twmoMultNTracksTPCOnlyUnfm80, float);
DECLARE_SOA_COLUMN(TwmoMultNTracksITSTPCUnfm80, twmoMultNTracksITSTPCUnfm80, float);
DECLARE_SOA_COLUMN(TwmoMultAllTracksTPCOnlyUnfm80, twmoMultAllTracksTPCOnlyUnfm80, float);

DECLARE_SOA_COLUMN(TwmoRobustT0V0PrimUnfm80, twmoRobustT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(TwmoRobustFDDT0V0PrimUnfm80, twmoRobustFDDT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(TwmoRobustNtrackDetUnfm80, twmoRobustNtrackDetUnfm80, float);
DECLARE_SOA_COLUMN(TwmoRobustMultExtraTableUnfm80, twmoRobustMultExtraTableUnfm80, float);

} // namespace trackmeanocc

// Tracks
// using Tracks = aod::Tracks;
// DECLARE_SOA_INDEX_TABLE_USER(TrackMeanOccs0, Tracks, "TRACKMEANOCCS0", o2::aod::trackmeanocc::TrackId);

DECLARE_SOA_TABLE(TmoTrackIds, "AOD", "TMOTRACKIDS", o2::aod::trackmeanocc::TrackId);

DECLARE_SOA_TABLE(TmoPrim, "AOD", "TMOPRIM", o2::soa::Index<>, // TrackMeanOccDet
                  o2::aod::trackmeanocc::TmoPrimUnfm80);

DECLARE_SOA_TABLE(TmoT0V0, "AOD", "TMOT0V0", o2::soa::Index<>, // TrackMeanOccDet
                  o2::aod::trackmeanocc::TmoFV0AUnfm80,
                  o2::aod::trackmeanocc::TmoFV0CUnfm80,
                  o2::aod::trackmeanocc::TmoFT0AUnfm80,
                  o2::aod::trackmeanocc::TmoFT0CUnfm80);

DECLARE_SOA_TABLE(TmoFDD, "AOD", "TMOFDD", o2::soa::Index<>, // TrackMeanOccDet
                  o2::aod::trackmeanocc::TmoFDDAUnfm80,
                  o2::aod::trackmeanocc::TmoFDDCUnfm80);

DECLARE_SOA_TABLE(TmoNTrackDet, "AOD", "TMONTRACKDET", o2::soa::Index<>, // TrackMeanOccNtrackDet
                  o2::aod::trackmeanocc::TmoNTrackITSUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackITSTPCCUnfm80);

DECLARE_SOA_TABLE(TmoMultExtra, "AOD", "TMOMULTEXTRA", o2::soa::Index<>, // TrackMeanOccMultExtra
                  o2::aod::trackmeanocc::TmoMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::TmoMultAllTracksTPCOnlyUnfm80);

DECLARE_SOA_TABLE(TmoRT0V0Prim, "AOD", "TMORT0V0PRIM", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TmoRobustT0V0PrimUnfm80);

DECLARE_SOA_TABLE(TmoRFDDT0V0Prim, "AOD", "TMORFDDT0V0PRIM", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TmoRobustFDDT0V0PrimUnfm80);

DECLARE_SOA_TABLE(TmoRNtrackDet, "AOD", "TMORNTRACKDET", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TmoRobustNtrackDetUnfm80);

DECLARE_SOA_TABLE(TmoRMultExtra, "AOD", "TMORMULTEXTRA", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TmoRobustMultExtraTableUnfm80);

DECLARE_SOA_TABLE(TwmoPrim, "AOD", "TWMOPRIM", o2::soa::Index<>, // WeightTrackMeanOcc
                  o2::aod::trackmeanocc::TwmoPrimUnfm80);

DECLARE_SOA_TABLE(TwmoT0V0, "AOD", "TWMOT0V0", o2::soa::Index<>, // WeightTrackMeanOccDet
                  o2::aod::trackmeanocc::TwmoFV0AUnfm80,
                  o2::aod::trackmeanocc::TwmoFV0CUnfm80,
                  o2::aod::trackmeanocc::TwmoFT0AUnfm80,
                  o2::aod::trackmeanocc::TwmoFT0CUnfm80);

DECLARE_SOA_TABLE(TwmoFDD, "AOD", "TWMOFDD", o2::soa::Index<>, // WeightTrackMeanOccDet
                  o2::aod::trackmeanocc::TwmoFDDAUnfm80,
                  o2::aod::trackmeanocc::TwmoFDDCUnfm80);

DECLARE_SOA_TABLE(TwmoNTrackDet, "AOD", "TWMONTRACKDET", o2::soa::Index<>, // WeightTrackMeanOccTrackMult
                  o2::aod::trackmeanocc::TwmoNTrackITSUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackITSTPCCUnfm80);

DECLARE_SOA_TABLE(TwmoMultExtra, "AOD", "TWMOMULTEXTRA", o2::soa::Index<>, // WeightTrackMeanOccMultExtra
                  o2::aod::trackmeanocc::TwmoMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoMultAllTracksTPCOnlyUnfm80);

DECLARE_SOA_TABLE(TwmoRT0V0Prim, "AOD", "TWMORT0V0PRIM", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TwmoRobustT0V0PrimUnfm80);

DECLARE_SOA_TABLE(TwmoRFDDT0V0Pri, "AOD", "TWMORFDDT0V0PRI", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TwmoRobustFDDT0V0PrimUnfm80);

DECLARE_SOA_TABLE(TwmoRNtrackDet, "AOD", "TWMORNTRACKDET", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TwmoRobustNtrackDetUnfm80);

DECLARE_SOA_TABLE(TwmoRMultExtra, "AOD", "TWMORMULTEXTRA", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TwmoRobustMultExtraTableUnfm80);

using Tmo = aod::TmoTrackIds::iterator;

using TrackQA = TracksQAVersion::iterator;

namespace trackmeanocc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Tmo, tmo, int64_t, TmoTrackIds, "");
DECLARE_SOA_INDEX_COLUMN_FULL(TrackQA, trackQA, int64_t, TracksQAVersion, "");
} // namespace trackmeanocc

DECLARE_SOA_TABLE(TrackToTracksQA, "AOD", "TRACKTOTRACKSQA", o2::aod::trackmeanocc::TrackQAId);
DECLARE_SOA_TABLE(TrackToTmo, "AOD", "TRACKTOTMO", o2::aod::trackmeanocc::TmoId);

DECLARE_SOA_TABLE(TrackQAToTmo, "AOD", "TRACKQATOTMO", o2::aod::trackmeanocc::TmoId);
DECLARE_SOA_TABLE(TmoToTrackQA, "AOD", "TMOTOTRACKQA", o2::aod::trackmeanocc::TrackQAId);

} // namespace o2::aod
#endif // COMMON_DATAMODEL_OCCUPANCYTABLES_H_
