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
namespace zdcdmoxygen
{
DECLARE_SOA_COLUMN(ZnaTdc, znaTdc, float);                 //! Tdc ZNA
DECLARE_SOA_COLUMN(ZnaAmpl, znaAmpl, float);               //! amplitude ZNA
DECLARE_SOA_COLUMN(ZnaPmc, znaPmc, float);                 //! ADC PmC ZNA
DECLARE_SOA_COLUMN(ZnaPm1, znaPm1, float);                 //! ADC Pm1 ZNA
DECLARE_SOA_COLUMN(ZnaPm2, znaPm2, float);                 //! ADC Pm2 ZNA
DECLARE_SOA_COLUMN(ZnaPm3, znaPm3, float);                 //! ADC Pm3 ZNA
DECLARE_SOA_COLUMN(ZnaPm4, znaPm4, float);                 //! ADC Pm4 ZNA
DECLARE_SOA_COLUMN(ZncTdc, zncTdc, float);                 //! Tdc ZNC
DECLARE_SOA_COLUMN(ZncAmpl, zncAmpl, float);               //! amplitude ZNC
DECLARE_SOA_COLUMN(ZncPmc, zncPmc, float);                 //! ADC PmC ZNC
DECLARE_SOA_COLUMN(ZncPm1, zncPm1, float);                 //! ADC Pm1 ZNC
DECLARE_SOA_COLUMN(ZncPm2, zncPm2, float);                 //! ADC Pm2 ZNC
DECLARE_SOA_COLUMN(ZncPm3, zncPm3, float);                 //! ADC Pm3 ZNC
DECLARE_SOA_COLUMN(ZncPm4, zncPm4, float);                 //! ADC Pm4 ZNC
DECLARE_SOA_COLUMN(ZpaTdc, zpaTdc, float);                 //! Tdc ZPA
DECLARE_SOA_COLUMN(ZpaAmpl, zpAmplc, float);               //! amplitude ZPA
DECLARE_SOA_COLUMN(ZpaPmc, zpaPmc, float);                 //! ADC PmC ZPA
DECLARE_SOA_COLUMN(ZpcTdc, zpcTdc, float);                 //! Tdc ZPC
DECLARE_SOA_COLUMN(ZpcAmpl, zpcAmpl, float);               //! amplitude ZPA
DECLARE_SOA_COLUMN(ZpcPmc, zpcPmc, float);                 //! ADC PmC ZPA
DECLARE_SOA_COLUMN(Zem1Tdc, zem1Tdc, float);               //! Tdc ZEM1
DECLARE_SOA_COLUMN(Zem1Ampl, zem1Ampl, float);              //! amplitude ZEM1
DECLARE_SOA_COLUMN(Zem2Tdc, zem2Tdc, float);               //! Tdc ZEM2
DECLARE_SOA_COLUMN(Zem2Ampl, zem2Ampl, float);             //! amplitude ZEM2
DECLARE_SOA_COLUMN(MultFt0a, multFt0a, float);             //! mult. FIT-A
DECLARE_SOA_COLUMN(MultFt0c, multFt0c, float);             //! mult. FIT-C
DECLARE_SOA_COLUMN(MultV0a, multV0a, float);               //! mult. V0-A
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);               //! Z vertex
DECLARE_SOA_COLUMN(CentralityFt0c, centralityFt0c, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFt0a, centralityFt0a, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFt0m, centralityFt0m, float); //! Centrality
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t); //! Selection Flags
} // namespace zdcdmoxygen

DECLARE_SOA_TABLE(ZDCdmOxygen, "AOD", "ZDCTABLELIGHTIONS",
                  zdcdmoxygen::ZnaTdc,
                  zdcdmoxygen::ZnaAmpl,
                  zdcdmoxygen::ZnaPmc,
                  zdcdmoxygen::ZnaPm1,
                  zdcdmoxygen::ZnaPm2,
                  zdcdmoxygen::ZnaPm3,
                  zdcdmoxygen::ZnaPm4,
                  zdcdmoxygen::ZncTdc,
                  zdcdmoxygen::ZncAmpl,
                  zdcdmoxygen::ZncPmc,
                  zdcdmoxygen::ZncPm1,
                  zdcdmoxygen::ZncPm2,
                  zdcdmoxygen::ZncPm3,
                  zdcdmoxygen::ZncPm4,
                  zdcdmoxygen::ZpaTdc,
                  zdcdmoxygen::ZpaAmpl,
                  zdcdmoxygen::ZpcTdc,
                  zdcdmoxygen::ZpcAmpl,
                  zdcdmoxygen::Zem1Tdc,
                  zdcdmoxygen::Zem1Ampl,
                  zdcdmoxygen::Zem2Tdc,
                  zdcdmoxygen::Zem2Ampl,
                  zdcdmoxygen::MultFt0a,
                  zdcdmoxygen::MultFt0c,
                  zdcdmoxygen::MultV0a,
                  zdcdmoxygen::VertexZ,
                  zdcdmoxygen::CentralityFt0c,
                  zdcdmoxygen::CentralityFt0a,
                  zdcdmoxygen::CentralityFt0m,
                  zdcdmoxygen::SelectionBits);
} // namespace o2::aod

#endif // PWGMM_DATAMODEL_ZDCDMOXYGEN_H_
