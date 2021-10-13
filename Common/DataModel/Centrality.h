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
DECLARE_SOA_COLUMN(CentEstV0M, centV0M, float);                           //! Centrality percentile estimated from V0C+V0A multiplicities
DECLARE_SOA_COLUMN(CentEstRun2SPDTracklets, centRun2SPDTracklets, float); //! Run2 centrality percentile estimated from SPD tracklets multiplicity
DECLARE_SOA_COLUMN(CentEstRun2CL0, centRun2CL0, float);                   //! Run2 centrality percentile estimated from CL0 multiplicity
DECLARE_SOA_COLUMN(CentEstRun2CL1, centRun2CL1, float);                   //! Run2 centrality percentile estimated from CL1 multiplicity
} // namespace cent
DECLARE_SOA_TABLE(CentV0Ms, "AOD", "CENTV0M", cent::CentEstV0M);                      //! V0M estimated centrality table
DECLARE_SOA_TABLE(CentRun2SPDs, "AOD", "CENTRUN2SPD", cent::CentEstRun2SPDTracklets); //! Run2 SPD tracklets estimated centrality table
DECLARE_SOA_TABLE(CentRun2CL0s, "AOD", "CENTRUN2CL0", cent::CentEstRun2CL0);          //! Run2 CL0 estimated centrality table
DECLARE_SOA_TABLE(CentRun2CL1s, "AOD", "CENTRUN2CL1", cent::CentEstRun2CL1);          //! Run2 CL1 estimated centrality table
using CentV0M = CentV0Ms::iterator;
using CentRun2SPD = CentRun2SPDs::iterator;
using CentRun2CL0 = CentRun2CL0s::iterator;
using CentRun2CL1 = CentRun2CL1s::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_CENTRALITY_H_
