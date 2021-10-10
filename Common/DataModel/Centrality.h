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
DECLARE_SOA_COLUMN(CentEstV0M, centV0M, float);                   //! Centrality percentile estimated from V0C+V0A multiplicities
DECLARE_SOA_COLUMN(CentEstSPDTracklets, centSPDTracklets, float); //! Centrality percentile estimated from SPD tracklets multiplicity
DECLARE_SOA_COLUMN(CentEstCL0, centCL0, float);                   //! Centrality percentile estimated from CL0 multiplicity
DECLARE_SOA_COLUMN(CentEstCL1, centCL1, float);                   //! Centrality percentile estimated from CL1 multiplicity
} // namespace cent
DECLARE_SOA_TABLE(CentV0Ms, "AOD", "CENTV0M", cent::CentEstV0M);          //! V0M estimated centrality table
DECLARE_SOA_TABLE(CentSPDs, "AOD", "CENTSPD", cent::CentEstSPDTracklets); //! SPD tracklets estimated centrality table
DECLARE_SOA_TABLE(CentCL0s, "AOD", "CENTCL0", cent::CentEstCL0);          //! CL0 estimated centrality table
DECLARE_SOA_TABLE(CentCL1s, "AOD", "CENTCL1", cent::CentEstCL1);          //! CL1 estimated centrality table
using CentV0M = CentV0Ms::iterator;
using CentSPD = CentSPDs::iterator;
using CentCL0 = CentCL0s::iterator;
using CentCL1 = CentCL1s::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_CENTRALITY_H_
