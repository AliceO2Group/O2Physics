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
DECLARE_SOA_COLUMN(CentRun2V0M, centRun2V0M, float);                   //! Run 2 cent. from V0C+V0A multiplicities
DECLARE_SOA_COLUMN(CentRun2V0A, centRun2V0A, float);                   //! Run 2 cent. from V0A multiplicities
DECLARE_SOA_COLUMN(CentRun2SPDTracklets, centRun2SPDTracklets, float); //! Run 2 cent. from SPD tracklets multiplicity
DECLARE_SOA_COLUMN(CentRun2SPDClusters, centRun2SPDClusters, float);   //! Run 2 cent. from SPD clusters multiplicity
DECLARE_SOA_COLUMN(CentRun2CL0, centRun2CL0, float);                   //! Run 2 cent. from CL0 multiplicity
DECLARE_SOA_COLUMN(CentRun2CL1, centRun2CL1, float);                   //! Run 2 cent. from CL1 multiplicity
DECLARE_SOA_COLUMN(CentRun2RefMult5, centRun2RefMult5, float);         //! Run 2 cent. from ref. mult. estimator, eta 0.5
DECLARE_SOA_COLUMN(CentRun2RefMult8, centRun2RefMult8, float);         //! Run 2 cent. from  ref. mult. estimator, eta 0.8

DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);                 //! Run 3 cent. from FV0A multiplicities
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);                 //! Run 3 cent. from FT0A+FT0C multiplicities
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);                 //! Run 3 cent. from FT0A multiplicity
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);                 //! Run 3 cent. from FT0C multiplicity
DECLARE_SOA_COLUMN(CentFT0CVariant1, centFT0CVariant1, float); //! Run 3 cent. from FT0C multiplicity
DECLARE_SOA_COLUMN(CentFDDM, centFDDM, float);                 //! Run 3 cent. from FDDA+FDDC multiplicity
DECLARE_SOA_COLUMN(CentNTPV, centNTPV, float);                 //! Run 3 cent. from the number of tracks contributing to the
DECLARE_SOA_COLUMN(CentNGlobal, centNGlobal, float);           //! Run 3 cent. from the number of tracks contributing to the PV
DECLARE_SOA_COLUMN(CentMFT, centMFT, float);                   //! Run 3 cent. from the number of tracks in the MFT
} // namespace cent

// Run 2 tables
DECLARE_SOA_TABLE(CentRun2V0Ms, "AOD", "CENTRUN2V0M", cent::CentRun2V0M);                //! Run 2 V0M centrality table
DECLARE_SOA_TABLE(CentRun2V0As, "AOD", "CENTRUN2V0A", cent::CentRun2V0A);                //! Run 2 V0A centrality table
DECLARE_SOA_TABLE(CentRun2SPDTrks, "AOD", "CENTRUN2SPDTRK", cent::CentRun2SPDTracklets); //! Run 2 SPD tracklets centrality table
DECLARE_SOA_TABLE(CentRun2SPDClss, "AOD", "CENTRUN2SPDCLS", cent::CentRun2SPDClusters);  //! Run 2 SPD clusters centrality table
DECLARE_SOA_TABLE(CentRun2CL0s, "AOD", "CENTRUN2CL0", cent::CentRun2CL0);                //! Run 2 CL0 centrality table
DECLARE_SOA_TABLE(CentRun2CL1s, "AOD", "CENTRUN2CL1", cent::CentRun2CL1);                //! Run 2 CL1 centrality table
DECLARE_SOA_TABLE(CentRun2RefMult5s, "AOD", "CENTRUN2REFMULT5", cent::CentRun2RefMult5); //! Run 2, ref mult |eta| < 0.5
DECLARE_SOA_TABLE(CentRun2RefMult8s, "AOD", "CENTRUN2REFMULT8", cent::CentRun2RefMult8); //! Run 2, ref mult |eta| < 0.8

// Run 3 tables
DECLARE_SOA_TABLE(CentFV0As, "AOD", "CENTFV0A", cent::CentFV0A);          //! Run 3 FV0A centrality table
DECLARE_SOA_TABLE(CentFT0Ms, "AOD", "CENTFT0M", cent::CentFT0M);          //! Run 3 FT0M centrality table
DECLARE_SOA_TABLE(CentFT0As, "AOD", "CENTFT0A", cent::CentFT0A);          //! Run 3 FT0A centrality table
DECLARE_SOA_TABLE(CentFT0Cs, "AOD", "CENTFT0C", cent::CentFT0C);          //! Run 3 FT0C centrality table
DECLARE_SOA_TABLE(CentFDDMs, "AOD", "CENTFDDM", cent::CentFDDM);          //! Run 3 FDDM centrality table
DECLARE_SOA_TABLE(CentNTPVs, "AOD", "CENTNTPV", cent::CentNTPV);          //! Run 3 NTPV centrality table
DECLARE_SOA_TABLE(CentNGlobals, "AOD", "CENTNGLOBAL", cent::CentNGlobal); //! Run 3 NGlobal centrality table
DECLARE_SOA_TABLE(CentMFTs, "AOD", "CENTMFT", cent::CentMFT);             //! Run 3 MFT tracks centrality table

// Run 3 variant tables
DECLARE_SOA_TABLE(CentFT0CVariant1s, "AOD", "CENTFT0Cvar1", cent::CentFT0CVariant1); //! Run 3 FT0C variant 1

using CentRun2V0M = CentRun2V0Ms::iterator;
using CentRun2V0A = CentRun2V0As::iterator;
using CentRun2SPDTrk = CentRun2SPDTrks::iterator;
using CentRun2SPDCls = CentRun2SPDClss::iterator;
using CentRun2CL0 = CentRun2CL0s::iterator;
using CentRun2CL1 = CentRun2CL1s::iterator;
using CentRun2RefMult5 = CentRun2RefMult5s::iterator;
using CentRun2RefMult8 = CentRun2RefMult8s::iterator;
using CentFV0A = CentFV0As::iterator;
using CentFT0M = CentFT0Ms::iterator;
using CentFT0A = CentFT0As::iterator;
using CentFT0C = CentFT0Cs::iterator;
using CentFDDM = CentFDDMs::iterator;
using CentNTPV = CentNTPVs::iterator;
using CentNGlobal = CentNGlobals::iterator;
using CentMFT = CentMFTs::iterator;

template <typename T>
concept HasRun2Centrality = requires(T&& t) {
  { t.centRun2V0M() };
  { t.centRun2CL0() };
  { t.centRun2CL1() };
};

template <typename T>
concept HasCentrality = requires(T&& t) {
  { t.centFV0A() };
  { t.centFT0M() };
  { t.centFT0A() };
  { t.centFT0C() };
  { t.centNTPV() };
};

} // namespace o2::aod

#endif // COMMON_DATAMODEL_CENTRALITY_H_
