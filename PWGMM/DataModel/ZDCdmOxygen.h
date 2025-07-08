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

/// \file ZDCInterCalib.h
/// \brief ZDC tower intercalibration task
/// \author Chiara Oppedisano <chiara.oppedisano@cern.ch>, INFN Torino

#ifndef COMMON_DATAMODEL_ZDCDMOXYGEN_H_
#define COMMON_DATAMODEL_ZDCDMOXYGEN_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace zdctable // o2-linter: disable=name/workflow-file
{
DECLARE_SOA_COLUMN(ZNAtdc, znaTDC, float);                 //! TDC ZNA           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNAampl, znaampl, float);               //! amplitude ZNA     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApmc, znapmc, float);                 //! ADC PMC ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm1, znaPM1, float);                 //! ADC PM1 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm2, znaPM2, float);                 //! ADC PM2 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm3, znaPM3, float);                 //! ADC PM3 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm4, znaPM4, float);                 //! ADC PM4 ZNA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCtdc, zncTDC, float);                 //! TDC ZNC           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCampl, zncampl, float);               //! amplitude ZNC     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpmc, zncpmc, float);                 //! ADC PMC ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm1, ZNCPM1, float);                 //! ADC PM1 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm2, ZNCPM2, float);                 //! ADC PM2 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm3, ZNCPM3, float);                 //! ADC PM3 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm4, ZNCPM4, float);                 //! ADC PM4 ZNC       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZPAtdc, zpaTDC, float);                 //! TDC ZPA           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZPAampl, zpamplc, float);               //! amplitude ZPA     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZPApmc, zpapmc, float);                 //! ADC PMC ZPA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZPCtdc, zpcTDC, float);                 //! TDC ZPC           // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZPCampl, zncampl, float);               //! amplitude ZPA     // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZPCpmc, zpcpmc, float);                 //! ADC PMC ZPA       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZEM1tdc, zem1TDC, float);               //! TDC ZEM1          // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZEM1ampl, zem1ampl, float);             //! amplitude ZEM1    // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZEM2tdc, zem2TDC, float);               //! TDC ZEM2          // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZEM2ampl, zem2ampl, float);             //! amplitude ZEM2    // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float);             //! mult. FIT-A       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);             //! mult. FIT-C       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(MultV0A, multV0A, float);               //! mult. V0-A       // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(CentralityFT0C, centralityFT0C, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFT0A, centralityFT0A, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFT0M, centralityFT0M, float); //! Centrality
DECLARE_SOA_COLUMN(CentralityFV0A, centralityFV0A, float); //! Centrality
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t); //! Selection Flags
} // namespace zdctable

DECLARE_SOA_TABLE(ZDCInterCalib, "AOD", "ZDCO", o2::soa::Index<>,
                  zdctable::ZNAtdc,
                  zdctable::ZNAampl,
                  zdctable::ZNApmc,
                  zdctable::ZNApm1,
                  zdctable::ZNApm2,
                  zdctable::ZNApm3,
                  zdctable::ZNApm4,
                  zdctable::ZNCtdc,
                  zdctable::ZNCampl,
                  zdctable::ZNCpmc,
                  zdctable::ZNCpm1,
                  zdctable::ZNCpm2,
                  zdctable::ZNCpm3,
                  zdctable::ZNCpm4,
                  zdctable::ZPAtdc,
                  zdctable::ZPAampl,
                  zdctable::ZPCtdc,
                  zdctable::ZPCampl,
                  zdctable::ZEM1tdc,
                  zdctable::ZEM1ampl
                    zdctable::ZEM2tdc,
                  zdctable::ZEM2ampl,
                  zdctable::MultFT0A,
                  zdctable::MultFT0C,
                  zdctable::MultV0A,
                  zdctable::CentralityFT0C,
                  zdctable::CentralityFT0A,
                  zdctable::CentralityFT0M,
                  zdctable::CentralityFV0A,
                  zdctable::SelectionBits);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_ZDCDMOXYGEN_H_
