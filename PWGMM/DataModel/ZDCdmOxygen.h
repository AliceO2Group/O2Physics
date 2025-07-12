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

/// \file ZDCdmOxygen.h
/// \brief ZDC tale for O-O Ne-Ne and p-O collisions
/// \author Chiara Oppedisano <chiara.oppedisano@cern.ch>, INFN Torino

#ifndef PWGMM_DATAMODEL_ZDCDMOXYGEN_H_
#define PWGMM_DATAMODEL_ZDCDMOXYGEN_H_

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Centrality.h"

namespace o2::aod
{
namespace zdcTableOO // o2-linter: disable=name/namespace
{
DECLARE_SOA_COLUMN(ZnaTdc, znaTDC, float);                 //! TDC ZNA           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZnaAmpl, znaAmpl, float);               //! amplitude ZNA     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZnaPmc, znaPMC, float);                 //! ADC PMC ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZnaPm1, znaPM1, float);                 //! ADC PM1 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZnaPm2, znaPM2, float);                 //! ADC PM2 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZnaPm3, znaPM3, float);                 //! ADC PM3 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZnaPm4, znaPM4, float);                 //! ADC PM4 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZncTdc, zncTDC, float);                 //! TDC ZNC           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZncAmpl, zncAmpl, float);               //! amplitude ZNC     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZncPmc, zncPMC, float);                 //! ADC PMC ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZncPm1, zncPM1, float);                 //! ADC PM1 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZncPm2, zncPM2, float);                 //! ADC PM2 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZncPm3, zncPM3, float);                 //! ADC PM3 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZncPm4, zncPM4, float);                 //! ADC PM4 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZpaTdc, zpaTDC, float);                 //! TDC ZPA           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZpaAmpl, zpAmplc, float);               //! amplitude ZPA     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZpaPmc, zpaPMC, float);                 //! ADC PMC ZPA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZpaTdc, zpcTDC, float);                 //! TDC ZPC           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZpaAmpl, zpcAmpl, float);               //! amplitude ZPA     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZpaPmc, zpcPMC, float);                 //! ADC PMC ZPA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(Zem1Tdc, zem1TDC, float);               //! TDC ZEM1          // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(Zem1Ampl, zemampl, float);              //! amplitude ZEM1    // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(Zem2Tdc, zem2TDC, float);               //! TDC ZEM2          // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(Zem2ampl, zem2Ampl, float);             //! amplitude ZEM2    // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(MultFt0A, multFT0A, float);             //! mult. FIT-A       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(MultFt0C, multFT0C, float);             //! mult. FIT-C       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(MultV0A, multV0A, float);               //! mult. V0-A        // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(Zvertex, zVertex, float);               //! Z vertex          // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(CentralityFt0C, centralityFT0C, float); //! Centrality        // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(CentralityFt0A, centralityFT0A, float); //! Centrality        // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(CentralityFt0M, centralityFT0M, float); //! Centrality        // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t); //! Selection Flags   // o2-linter: disable=name/o2-column
} // namespace zdctable

DECLARE_SOA_TABLE(ZdcTable, "AOD", "ZdcTeble",
                  zdcTableOO::ZNAtdc,
                  zdcTableOO::ZNAampl,
                  zdcTableOO::ZNApmc,
                  zdcTableOO::ZNApm1,
                  zdcTableOO::ZNApm2,
                  zdcTableOO::ZNApm3,
                  zdcTableOO::ZNApm4,
                  zdcTableOO::ZNCtdc,
                  zdcTableOO::ZNCampl,
                  zdcTableOO::ZNCpmc,
                  zdcTableOO::ZNCpm1,
                  zdcTableOO::ZNCpm2,
                  zdcTableOO::ZNCpm3,
                  zdcTableOO::ZNCpm4,
                  zdcTableOO::ZPAtdc,
                  zdcTableOO::ZPAampl,
                  zdcTableOO::ZPCtdc,
                  zdcTableOO::ZPCampl,
                  zdcTableOO::ZEM1tdc,
                  zdcTableOO::ZEM1ampl,
                  zdcTableOO::ZEM2tdc,
                  zdcTableOO::ZEM2ampl,
                  zdcTableOO::MultFT0A,
                  zdcTableOO::MultFT0C,
                  zdcTableOO::MultV0A,
                  zdcTableOO::Zvertex,
                  zdcTableOO::CentralityFT0C,
                  zdcTableOO::CentralityFT0A,
                  zdcTableOO::CentralityFT0M,
                  zdcTableOO::SelectionBits);
} // namespace o2::aod

#endif // PWGMM_DATAMODEL_ZDCDMOXYGEN_H_
