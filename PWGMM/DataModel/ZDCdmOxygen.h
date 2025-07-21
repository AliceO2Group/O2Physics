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

#include "Common/DataModel/Centrality.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace zdcTableLiIo
{
DECLARE_SOA_COLUMN(ZnaTdc, znaTdc, float);                 //! Tdc ZNA
DECLARE_SOA_COLUMN(ZnaAmpl, znaAmpl, float);               //! amplitude ZNA
DECLARE_SOA_COLUMN(ZnaPmc, znaPmC, float);                 //! ADC PmC ZNA
DECLARE_SOA_COLUMN(ZnaPm1, znaPm1, float);                 //! ADC Pm1 ZNA
DECLARE_SOA_COLUMN(ZnaPm2, znaPm2, float);                 //! ADC Pm2 ZNA
DECLARE_SOA_COLUMN(ZnaPm3, znaPm3, float);                 //! ADC Pm3 ZNA
DECLARE_SOA_COLUMN(ZnaPm4, znaPm4, float);                 //! ADC Pm4 ZNA
DECLARE_SOA_COLUMN(ZncTdc, zncTdc, float);                 //! Tdc ZNC
DECLARE_SOA_COLUMN(ZncAmpl, zncAmpl, float);               //! amplitude ZNC
DECLARE_SOA_COLUMN(ZncPmc, zncPmC, float);                 //! ADC PmC ZNC
DECLARE_SOA_COLUMN(ZncPm1, zncPm1, float);                 //! ADC Pm1 ZNC
DECLARE_SOA_COLUMN(ZncPm2, zncPm2, float);                 //! ADC Pm2 ZNC
DECLARE_SOA_COLUMN(ZncPm3, zncPm3, float);                 //! ADC Pm3 ZNC
DECLARE_SOA_COLUMN(ZncPm4, zncPm4, float);                 //! ADC Pm4 ZNC
DECLARE_SOA_COLUMN(ZpaTdc, zpaTdc, float);                 //! Tdc ZPA
DECLARE_SOA_COLUMN(ZpaAmpl, zpAmplc, float);               //! amplitude ZPA
DECLARE_SOA_COLUMN(ZpaPmc, zpaPmC, float);                 //! ADC PmC ZPA
DECLARE_SOA_COLUMN(ZpcTdc, zpcTdc, float);                 //! Tdc ZPC
DECLARE_SOA_COLUMN(ZpcAmpl, zpcAmpl, float);               //! amplitude ZPA
DECLARE_SOA_COLUMN(ZpcPmc, zpcPmC, float);                 //! ADC PmC ZPA
DECLARE_SOA_COLUMN(Zem1Tdc, zem1Tdc, float);               //! Tdc ZEM1
DECLARE_SOA_COLUMN(Zem1Ampl, zemampl, float);              //! amplitude ZEM1
DECLARE_SOA_COLUMN(Zem2Tdc, zem2Tdc, float);               //! Tdc ZEM2
DECLARE_SOA_COLUMN(Zem2Ampl, zem2Ampl, float);             //! amplitude ZEM2
DECLARE_SOA_COLUMN(MultFt0A, multFT0A, float);             //! mult. FIT-A
DECLARE_SOA_COLUMN(MultFt0C, multFT0C, float);             //! mult. FIT-C
DECLARE_SOA_COLUMN(MultV0A, multV0A, float);               //! mult. V0-A
DECLARE_SOA_COLUMN(Zvertex, zVertex, float);               //! Z vertex
DECLARE_SOA_COLUMN(CentralityFt0C, centralityFT0C, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFt0A, centralityFT0A, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFt0M, centralityFT0M, float); //! Centrality
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t); //! Selection Flags
} // namespace zdcTableLiIo

DECLARE_SOA_TABLE(ZdcTable, "AOD", "ZdcTeble",
                  zdcTableLiIo::ZnaTdc,
                  zdcTableLiIo::ZnaAmpl,
                  zdcTableLiIo::ZnaPmc,
                  zdcTableLiIo::ZnaPm1,
                  zdcTableLiIo::ZnaPm2,
                  zdcTableLiIo::ZnaPm3,
                  zdcTableLiIo::ZnaPm4,
                  zdcTableLiIo::ZncTdc,
                  zdcTableLiIo::ZncAmpl,
                  zdcTableLiIo::ZncPmc,
                  zdcTableLiIo::ZncPm1,
                  zdcTableLiIo::ZncPm2,
                  zdcTableLiIo::ZncPm3,
                  zdcTableLiIo::ZncPm4,
                  zdcTableLiIo::ZpaTdc,
                  zdcTableLiIo::ZpaAmpl,
                  zdcTableLiIo::ZpcTdc,
                  zdcTableLiIo::ZpcAmpl,
                  zdcTableLiIo::Zem1Tdc,
                  zdcTableLiIo::Zem1Ampl,
                  zdcTableLiIo::Zem2Tdc,
                  zdcTableLiIo::Zem2Ampl,
                  zdcTableLiIo::MultFt0A,
                  zdcTableLiIo::MultFt0C,
                  zdcTableLiIo::MultV0A,
                  zdcTableLiIo::Zvertex,
                  zdcTableLiIo::CentralityFt0C,
                  zdcTableLiIo::CentralityFt0A,
                  zdcTableLiIo::CentralityFt0M,
                  zdcTableLiIo::SelectionBits);
} // namespace o2::aod

#endif // PWGMM_DATAMODEL_ZDCDMOXYGEN_H_
