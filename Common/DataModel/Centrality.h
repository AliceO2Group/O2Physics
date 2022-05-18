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
#ifndef O2_ANALYSIS_CENTRALITY_H_
#define O2_ANALYSIS_CENTRALITY_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace cent
{
DECLARE_SOA_COLUMN(CentRun2V0M, centRun2V0M, float);                   //! Centrality percentile estimated from V0C+V0A multiplicities
DECLARE_SOA_COLUMN(CentRun2SPDTracklets, centRun2SPDTracklets, float); //! Run2 centrality percentile estimated from SPD tracklets multiplicity
DECLARE_SOA_COLUMN(CentRun2SPDClusters, centRun2SPDClusters, float);   //! Run2 centrality percentile estimated from SPD clusters multiplicity
DECLARE_SOA_COLUMN(CentRun2CL0, centRun2CL0, float);                   //! Run2 centrality percentile estimated from CL0 multiplicity
DECLARE_SOA_COLUMN(CentRun2CL1, centRun2CL1, float);                   //! Run2 centrality percentile estimated from CL1 multiplicity
} // namespace cent
DECLARE_SOA_TABLE(CentRun2V0Ms, "AOD", "CENTRUN2V0M", cent::CentRun2V0M);                //! V0M estimated centrality table
DECLARE_SOA_TABLE(CentRun2SPDTrks, "AOD", "CENTRUN2SPDTRK", cent::CentRun2SPDTracklets); //! Run2 SPD tracklets estimated centrality table
DECLARE_SOA_TABLE(CentRun2SPDClss, "AOD", "CENTRUN2SPDCLS", cent::CentRun2SPDClusters);  //! Run2 SPD clusters estimated centrality table
DECLARE_SOA_TABLE(CentRun2CL0s, "AOD", "CENTRUN2CL0", cent::CentRun2CL0);                //! Run2 CL0 estimated centrality table
DECLARE_SOA_TABLE(CentRun2CL1s, "AOD", "CENTRUN2CL1", cent::CentRun2CL1);                //! Run2 CL1 estimated centrality table
using CentRun2V0M = CentRun2V0Ms::iterator;
using CentRun2SPDTrk = CentRun2SPDTrks::iterator;
using CentRun2SPDCls = CentRun2SPDClss::iterator;
using CentRun2CL0 = CentRun2CL0s::iterator;
using CentRun2CL1 = CentRun2CL1s::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_CENTRALITY_H_
