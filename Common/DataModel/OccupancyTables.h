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
DECLARE_SOA_COLUMN(BCsInTFlist, bcsInTFlist, std::vector<int64_t>);

DECLARE_SOA_COLUMN(Occ_Prim_Unfm_80, occ_Prim_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_FV0A_Unfm_80, occ_FV0A_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_FV0C_Unfm_80, occ_FV0C_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_FT0A_Unfm_80, occ_FT0A_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_FT0C_Unfm_80, occ_FT0C_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_FDDA_Unfm_80, occ_FDDA_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_FDDC_Unfm_80, occ_FDDC_Unfm_80, std::vector<float>);

DECLARE_SOA_COLUMN(Occ_NTrack_PVC_Unfm_80, occ_NTrack_PVC_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrack_ITS_Unfm_80, occ_NTrack_ITS_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrack_TPC_Unfm_80, occ_NTrack_TPC_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrack_TRD_Unfm_80, occ_NTrack_TRD_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrack_TOF_Unfm_80, occ_NTrack_TOF_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrackSize_Unfm_80, occ_NTrackSize_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrackTPC_A_Unfm_80, occ_NTrackTPC_A_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrackTPC_C_Unfm_80, occ_NTrackTPC_C_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrackITS_TPC_Unfm_80, occ_NTrackITS_TPC_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrackITS_TPC_A_Unfm_80, occ_NTrackITS_TPC_A_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(Occ_NTrackITS_TPC_C_Unfm_80, occ_NTrackITS_TPC_C_Unfm_80, std::vector<float>);

DECLARE_SOA_COLUMN(OccRobust_T0V0Prim_Unfm_80, occRobust_T0V0Prim_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(OccRobust_FDDT0V0Prim_Unfm_80, occRobust_FDDT0V0Prim_Unfm_80, std::vector<float>);
DECLARE_SOA_COLUMN(OccRobust_NtrackDet_Unfm_80, occRobust_NtrackDet_Unfm_80, std::vector<float>);

DECLARE_SOA_COLUMN(MeanOcc_Prim_Unfm_80, meanOcc_Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FV0A_Unfm_80, meanOcc_FV0A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FV0C_Unfm_80, meanOcc_FV0C_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FT0A_Unfm_80, meanOcc_FT0A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FT0C_Unfm_80, meanOcc_FT0C_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FDDA_Unfm_80, meanOcc_FDDA_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FDDC_Unfm_80, meanOcc_FDDC_Unfm_80, float);

DECLARE_SOA_COLUMN(MeanOcc_NTrack_PVC_Unfm_80, meanOcc_NTrack_PVC_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_ITS_Unfm_80, meanOcc_NTrack_ITS_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_TPC_Unfm_80, meanOcc_NTrack_TPC_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_TRD_Unfm_80, meanOcc_NTrack_TRD_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_TOF_Unfm_80, meanOcc_NTrack_TOF_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackSize_Unfm_80, meanOcc_NTrackSize_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackTPC_A_Unfm_80, meanOcc_NTrackTPC_A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackTPC_C_Unfm_80, meanOcc_NTrackTPC_C_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackITS_TPC_Unfm_80, meanOcc_NTrackITS_TPC_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackITS_TPC_A_Unfm_80, meanOcc_NTrackITS_TPC_A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackITS_TPC_C_Unfm_80, meanOcc_NTrackITS_TPC_C_Unfm_80, float);

DECLARE_SOA_COLUMN(MeanOccRobust_T0V0Prim_Unfm_80, meanOccRobust_T0V0Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOccRobust_FDDT0V0Prim_Unfm_80, meanOccRobust_FDDT0V0Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOccRobust_NtrackDet_Unfm_80, meanOccRobust_NtrackDet_Unfm_80, float);
} // namespace occp

DECLARE_SOA_TABLE(Occs, "AOD", "OCCS", o2::soa::Index<>, o2::aod::occp::TfId, o2::aod::occp::BCsInTFlist,
                  o2::aod::occp::Occ_Prim_Unfm_80,
                  o2::aod::occp::Occ_FV0A_Unfm_80,
                  o2::aod::occp::Occ_FV0C_Unfm_80,
                  o2::aod::occp::Occ_FT0A_Unfm_80,
                  o2::aod::occp::Occ_FT0C_Unfm_80,
                  o2::aod::occp::Occ_FDDA_Unfm_80,
                  o2::aod::occp::Occ_FDDC_Unfm_80,
                  o2::aod::occp::Occ_NTrack_PVC_Unfm_80,
                  o2::aod::occp::Occ_NTrack_ITS_Unfm_80,
                  o2::aod::occp::Occ_NTrack_TPC_Unfm_80,
                  o2::aod::occp::Occ_NTrack_TRD_Unfm_80,
                  o2::aod::occp::Occ_NTrack_TOF_Unfm_80,
                  o2::aod::occp::Occ_NTrackSize_Unfm_80,
                  o2::aod::occp::Occ_NTrackTPC_A_Unfm_80,
                  o2::aod::occp::Occ_NTrackTPC_C_Unfm_80,
                  o2::aod::occp::Occ_NTrackITS_TPC_Unfm_80,
                  o2::aod::occp::Occ_NTrackITS_TPC_A_Unfm_80,
                  o2::aod::occp::Occ_NTrackITS_TPC_C_Unfm_80,
                  o2::aod::occp::OccRobust_T0V0Prim_Unfm_80,
                  o2::aod::occp::OccRobust_FDDT0V0Prim_Unfm_80,
                  o2::aod::occp::OccRobust_NtrackDet_Unfm_80,
                  o2::aod::occp::MeanOcc_Prim_Unfm_80,
                  o2::aod::occp::MeanOcc_FV0A_Unfm_80,
                  o2::aod::occp::MeanOcc_FV0C_Unfm_80,
                  o2::aod::occp::MeanOcc_FT0A_Unfm_80,
                  o2::aod::occp::MeanOcc_FT0C_Unfm_80,
                  o2::aod::occp::MeanOcc_FDDA_Unfm_80,
                  o2::aod::occp::MeanOcc_FDDC_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrack_PVC_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrack_ITS_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrack_TPC_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrack_TRD_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrack_TOF_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrackSize_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrackTPC_A_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrackTPC_C_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrackITS_TPC_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrackITS_TPC_A_Unfm_80,
                  o2::aod::occp::MeanOcc_NTrackITS_TPC_C_Unfm_80,
                  o2::aod::occp::MeanOccRobust_T0V0Prim_Unfm_80,
                  o2::aod::occp::MeanOccRobust_FDDT0V0Prim_Unfm_80,
                  o2::aod::occp::MeanOccRobust_NtrackDet_Unfm_80);

using Occ = Occs::iterator;
namespace occidx
{
DECLARE_SOA_INDEX_COLUMN(BC, bc);
DECLARE_SOA_INDEX_COLUMN(Occ, occ); // Iterator is passed here
DECLARE_SOA_COLUMN(TfId, tfId, int);
DECLARE_SOA_COLUMN(BCinTF, bcInTF, int);
} // namespace occidx

DECLARE_SOA_TABLE(OccIndexTable, "AOD", "OCCINDEXTABLE", o2::soa::Index<>,
                  o2::aod::occidx::BCId,
                  o2::aod::occidx::OccId,
                  o2::aod::occidx::TfId,
                  o2::aod::occidx::BCinTF);

namespace trackmeanocc
{
DECLARE_SOA_INDEX_COLUMN(Track, track);

DECLARE_SOA_COLUMN(MeanOcc_Prim_Unfm_80, meanOcc_Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FV0A_Unfm_80, meanOcc_FV0A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FV0C_Unfm_80, meanOcc_FV0C_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FT0A_Unfm_80, meanOcc_FT0A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FT0C_Unfm_80, meanOcc_FT0C_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FDDA_Unfm_80, meanOcc_FDDA_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_FDDC_Unfm_80, meanOcc_FDDC_Unfm_80, float);

DECLARE_SOA_COLUMN(MeanOcc_NTrack_PVC_Unfm_80, meanOcc_NTrack_PVC_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_ITS_Unfm_80, meanOcc_NTrack_ITS_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_TPC_Unfm_80, meanOcc_NTrack_TPC_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_TRD_Unfm_80, meanOcc_NTrack_TRD_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrack_TOF_Unfm_80, meanOcc_NTrack_TOF_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackSize_Unfm_80, meanOcc_NTrackSize_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackTPC_A_Unfm_80, meanOcc_NTrackTPC_A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackTPC_C_Unfm_80, meanOcc_NTrackTPC_C_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackITS_TPC_Unfm_80, meanOcc_NTrackITS_TPC_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackITS_TPC_A_Unfm_80, meanOcc_NTrackITS_TPC_A_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOcc_NTrackITS_TPC_C_Unfm_80, meanOcc_NTrackITS_TPC_C_Unfm_80, float);

DECLARE_SOA_COLUMN(MeanOccRobust_T0V0Prim_Unfm_80, meanOccRobust_T0V0Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOccRobust_FDDT0V0Prim_Unfm_80, meanOccRobust_FDDT0V0Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(MeanOccRobust_NtrackDet_Unfm_80, meanOccRobust_NtrackDet_Unfm_80, float);

DECLARE_SOA_COLUMN(WeightMeanOcc_Prim_Unfm_80, weightMeanOcc_Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_FV0A_Unfm_80, weightMeanOcc_FV0A_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_FV0C_Unfm_80, weightMeanOcc_FV0C_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_FT0A_Unfm_80, weightMeanOcc_FT0A_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_FT0C_Unfm_80, weightMeanOcc_FT0C_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_FDDA_Unfm_80, weightMeanOcc_FDDA_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_FDDC_Unfm_80, weightMeanOcc_FDDC_Unfm_80, float);

DECLARE_SOA_COLUMN(WeightMeanOcc_NTrack_PVC_Unfm_80, weightMeanOcc_NTrack_PVC_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrack_ITS_Unfm_80, weightMeanOcc_NTrack_ITS_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrack_TPC_Unfm_80, weightMeanOcc_NTrack_TPC_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrack_TRD_Unfm_80, weightMeanOcc_NTrack_TRD_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrack_TOF_Unfm_80, weightMeanOcc_NTrack_TOF_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrackSize_Unfm_80, weightMeanOcc_NTrackSize_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrackTPC_A_Unfm_80, weightMeanOcc_NTrackTPC_A_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrackTPC_C_Unfm_80, weightMeanOcc_NTrackTPC_C_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrackITS_TPC_Unfm_80, weightMeanOcc_NTrackITS_TPC_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrackITS_TPC_A_Unfm_80, weightMeanOcc_NTrackITS_TPC_A_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOcc_NTrackITS_TPC_C_Unfm_80, weightMeanOcc_NTrackITS_TPC_C_Unfm_80, float);

DECLARE_SOA_COLUMN(WeightMeanOccRobust_T0V0Prim_Unfm_80, weightMeanOccRobust_T0V0Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOccRobust_FDDT0V0Prim_Unfm_80, weightMeanOccRobust_FDDT0V0Prim_Unfm_80, float);
DECLARE_SOA_COLUMN(WeightMeanOccRobust_NtrackDet_Unfm_80, weightMeanOccRobust_NtrackDet_Unfm_80, float);

} // namespace trackmeanocc

DECLARE_SOA_TABLE(TrackMeanOccs, "AOD", "TRACKMEANOCCS", o2::soa::Index<>,
                  o2::aod::trackmeanocc::TrackId,
                  o2::aod::trackmeanocc::MeanOcc_Prim_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_FV0A_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_FV0C_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_FT0A_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_FT0C_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_FDDA_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_FDDC_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrack_PVC_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrack_ITS_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrack_TPC_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrack_TRD_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrack_TOF_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrackSize_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrackTPC_A_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrackTPC_C_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrackITS_TPC_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrackITS_TPC_A_Unfm_80,
                  o2::aod::trackmeanocc::MeanOcc_NTrackITS_TPC_C_Unfm_80,
                  o2::aod::trackmeanocc::MeanOccRobust_T0V0Prim_Unfm_80,
                  o2::aod::trackmeanocc::MeanOccRobust_FDDT0V0Prim_Unfm_80,
                  o2::aod::trackmeanocc::MeanOccRobust_NtrackDet_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_Prim_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_FV0A_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_FV0C_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_FT0A_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_FT0C_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_FDDA_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_FDDC_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrack_PVC_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrack_ITS_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrack_TPC_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrack_TRD_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrack_TOF_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrackSize_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrackTPC_A_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrackTPC_C_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrackITS_TPC_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrackITS_TPC_A_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOcc_NTrackITS_TPC_C_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOccRobust_T0V0Prim_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOccRobust_FDDT0V0Prim_Unfm_80,
                  o2::aod::trackmeanocc::WeightMeanOccRobust_NtrackDet_Unfm_80);
} // namespace o2::aod
#endif // COMMON_DATAMODEL_OCCUPANCYTABLES_H_
