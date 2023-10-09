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
#ifndef COMMON_DATAMODEL_CENTRALITY_H_
#define COMMON_DATAMODEL_CENTRALITY_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace cent
{
DECLARE_SOA_COLUMN(CentRun2V0M, centRun2V0M, float);                   //! Run2 Centrality percentile estimated from V0C+V0A multiplicities
DECLARE_SOA_COLUMN(CentRun2V0A, centRun2V0A, float);                   //! Run2 Centrality percentile estimated from V0A multiplicities
DECLARE_SOA_COLUMN(CentRun2SPDTracklets, centRun2SPDTracklets, float); //! Run2 centrality percentile estimated from SPD tracklets multiplicity
DECLARE_SOA_COLUMN(CentRun2SPDClusters, centRun2SPDClusters, float);   //! Run2 centrality percentile estimated from SPD clusters multiplicity
DECLARE_SOA_COLUMN(CentRun2CL0, centRun2CL0, float);                   //! Run2 centrality percentile estimated from CL0 multiplicity
DECLARE_SOA_COLUMN(CentRun2CL1, centRun2CL1, float);                   //! Run2 centrality percentile estimated from CL1 multiplicity
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);                         //! Run3 Centrality percentile estimated from FV0A multiplicities
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);                         //! Run3 centrality percentile estimated from FT0A+FT0C multiplicities
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);                         //! Run3 centrality percentile estimated from FT0A multiplicity
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);                         //! Run3 centrality percentile estimated from FT0C multiplicity
DECLARE_SOA_COLUMN(CentFDDM, centFDDM, float);                         //! Run3 centrality percentile estimated from FDDA+FDDC multiplicity
DECLARE_SOA_COLUMN(CentNTPV, centNTPV, float);                         //! Run3 centrality percentile estimated from the number of tracks contributing to the PV
} // namespace cent
DECLARE_SOA_TABLE(CentRun2V0Ms, "AOD", "CENTRUN2V0M", cent::CentRun2V0M);                //! Run2 V0M estimated centrality table
DECLARE_SOA_TABLE(CentRun2V0As, "AOD", "CENTRUN2V0A", cent::CentRun2V0A);                //! Run2 V0A estimated centrality table
DECLARE_SOA_TABLE(CentRun2SPDTrks, "AOD", "CENTRUN2SPDTRK", cent::CentRun2SPDTracklets); //! Run2 SPD tracklets estimated centrality table
DECLARE_SOA_TABLE(CentRun2SPDClss, "AOD", "CENTRUN2SPDCLS", cent::CentRun2SPDClusters);  //! Run2 SPD clusters estimated centrality table
DECLARE_SOA_TABLE(CentRun2CL0s, "AOD", "CENTRUN2CL0", cent::CentRun2CL0);                //! Run2 CL0 estimated centrality table
DECLARE_SOA_TABLE(CentRun2CL1s, "AOD", "CENTRUN2CL1", cent::CentRun2CL1);                //! Run2 CL1 estimated centrality table
DECLARE_SOA_TABLE(CentFV0As, "AOD", "CENTFV0A", cent::CentFV0A);                         //! Run3 FV0A estimated centrality table
DECLARE_SOA_TABLE(CentFT0Ms, "AOD", "CENTFT0M", cent::CentFT0M);                         //! Run3 FT0M estimated centrality table
DECLARE_SOA_TABLE(CentFT0As, "AOD", "CENTFT0A", cent::CentFT0A);                         //! Run3 FT0A estimated centrality table
DECLARE_SOA_TABLE(CentFT0Cs, "AOD", "CENTFT0C", cent::CentFT0C);                         //! Run3 FT0C estimated centrality table
DECLARE_SOA_TABLE(CentFDDMs, "AOD", "CENTFDDM", cent::CentFDDM);                         //! Run3 FDDM estimated centrality table
DECLARE_SOA_TABLE(CentNTPVs, "AOD", "CENTNTPV", cent::CentNTPV);                         //! Run3 NTPV estimated centrality table
using CentRun2V0M = CentRun2V0Ms::iterator;
using CentRun2V0A = CentRun2V0As::iterator;
using CentRun2SPDTrk = CentRun2SPDTrks::iterator;
using CentRun2SPDCls = CentRun2SPDClss::iterator;
using CentRun2CL0 = CentRun2CL0s::iterator;
using CentRun2CL1 = CentRun2CL1s::iterator;
using CentFV0A = CentFV0As::iterator;
using CentFT0M = CentFT0Ms::iterator;
using CentFT0A = CentFT0As::iterator;
using CentFT0C = CentFT0Cs::iterator;
using CentFDDM = CentFDDMs::iterator;
using CentNTPV = CentNTPVs::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_CENTRALITY_H_
