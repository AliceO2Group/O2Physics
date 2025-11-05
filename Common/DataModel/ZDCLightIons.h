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

/// \file ZDCLightIons.h
/// \brief ZDC data model for O-O Ne-Ne and p-O collisions
/// \author Chiara Oppedisano <chiara.oppedisano@cern.ch>

#ifndef PWGMM_DATAMODEL_ZDCLIGHTIONS_H_
#define PWGMM_DATAMODEL_ZDCLIGHTIONS_H_

#include "Common/DataModel/Centrality.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace zdclightions
{
DECLARE_SOA_COLUMN(ZnaTdc, znaTdc, float);                 //! TDC ZNA
DECLARE_SOA_COLUMN(ZnaAmpl, znaAmpl, float);               //! amplitude ZNA
DECLARE_SOA_COLUMN(ZnaPmc, znaPmc, float);                 //! ADC PmC ZNA
DECLARE_SOA_COLUMN(ZnaPm1, znaPm1, float);                 //! ADC PmQ1 ZNA
DECLARE_SOA_COLUMN(ZnaPm2, znaPm2, float);                 //! ADC PmQ2 ZNA
DECLARE_SOA_COLUMN(ZnaPm3, znaPm3, float);                 //! ADC PmQ3 ZNA
DECLARE_SOA_COLUMN(ZnaPm4, znaPm4, float);                 //! ADC PmQ4 ZNA
DECLARE_SOA_COLUMN(ZncTdc, zncTdc, float);                 //! TDC ZNC
DECLARE_SOA_COLUMN(ZncAmpl, zncAmpl, float);               //! amplitude ZNC
DECLARE_SOA_COLUMN(ZncPmc, zncPmc, float);                 //! ADC PmC ZNC
DECLARE_SOA_COLUMN(ZncPm1, zncPm1, float);                 //! ADC PmQ1 ZNC
DECLARE_SOA_COLUMN(ZncPm2, zncPm2, float);                 //! ADC PmQ2 ZNC
DECLARE_SOA_COLUMN(ZncPm3, zncPm3, float);                 //! ADC PmQ3 ZNC
DECLARE_SOA_COLUMN(ZncPm4, zncPm4, float);                 //! ADC PmQ4 ZNC
DECLARE_SOA_COLUMN(ZpaTdc, zpaTdc, float);                 //! TDC ZPA
DECLARE_SOA_COLUMN(ZpaAmpl, zpaAmpl, float);               //! amplitude ZPA
DECLARE_SOA_COLUMN(ZpaPmc, zpaPmc, float);                 //! ADC PmC ZPA
DECLARE_SOA_COLUMN(ZpcTdc, zpcTdc, float);                 //! TDC ZPC
DECLARE_SOA_COLUMN(ZpcAmpl, zpcAmpl, float);               //! amplitude ZPA
DECLARE_SOA_COLUMN(ZpcPmc, zpcPmc, float);                 //! ADC PmC ZPA
DECLARE_SOA_COLUMN(Zem1Tdc, zem1Tdc, float);               //! TDC ZEM1
DECLARE_SOA_COLUMN(Zem1Ampl, zem1Ampl, float);             //! amplitude ZEM1
DECLARE_SOA_COLUMN(Zem2Tdc, zem2Tdc, float);               //! TDC ZEM2
DECLARE_SOA_COLUMN(Zem2Ampl, zem2Ampl, float);             //! amplitude ZEM2
DECLARE_SOA_COLUMN(MultFt0a, multFt0a, float);             //! mult. FIT-A
DECLARE_SOA_COLUMN(MultFt0c, multFt0c, float);             //! mult. FIT-C
DECLARE_SOA_COLUMN(MultV0a, multV0a, float);               //! mult. V0-A
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);               //! Z vertex
DECLARE_SOA_COLUMN(CentralityFt0c, centralityFt0c, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFt0a, centralityFt0a, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFt0m, centralityFt0m, float); //! Centrality
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);        //! Timestamp
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t); //! Selection Flags
} // namespace zdclightions

DECLARE_SOA_TABLE(ZDCLightIons, "AOD", "ZDCTABLELIGHTIONS",
                  zdclightions::ZnaTdc,
                  zdclightions::ZnaAmpl,
                  zdclightions::ZnaPmc,
                  zdclightions::ZnaPm1,
                  zdclightions::ZnaPm2,
                  zdclightions::ZnaPm3,
                  zdclightions::ZnaPm4,
                  zdclightions::ZncTdc,
                  zdclightions::ZncAmpl,
                  zdclightions::ZncPmc,
                  zdclightions::ZncPm1,
                  zdclightions::ZncPm2,
                  zdclightions::ZncPm3,
                  zdclightions::ZncPm4,
                  zdclightions::ZpaTdc,
                  zdclightions::ZpaAmpl,
                  zdclightions::ZpcTdc,
                  zdclightions::ZpcAmpl,
                  zdclightions::Zem1Tdc,
                  zdclightions::Zem1Ampl,
                  zdclightions::Zem2Tdc,
                  zdclightions::Zem2Ampl,
                  zdclightions::MultFt0a,
                  zdclightions::MultFt0c,
                  zdclightions::MultV0a,
                  zdclightions::VertexZ,
                  zdclightions::CentralityFt0c,
                  zdclightions::CentralityFt0a,
                  zdclightions::CentralityFt0m,
                  zdclightions::Timestamp,
                  zdclightions::SelectionBits);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_ZDCLIGHTIONS_H_
