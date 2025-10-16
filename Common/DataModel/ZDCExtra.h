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

/// \file ZDCExtra.h
/// \brief ZDC extra table producer
/// \author Chiara Oppedisano <chiara.oppedisano@cern.ch>, INFN Torino
/// \author Uliana Dmitrieva <uliana.dmitrieva@cern.ch>, INFN Torino

#ifndef COMMON_DATAMODEL_ZDCEXTRA_H_
#define COMMON_DATAMODEL_ZDCEXTRA_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace zdcextra
{
DECLARE_SOA_COLUMN(ZNApmc, commonPMZNA, float);            //! PMC ZNA         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm1, ZNAPM1, float);                 //! PM1 ZNA         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm2, ZNAPM2, float);                 //! PM2 ZNA         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm3, ZNAPM3, float);                 //! PM3 ZNA         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm4, ZNAPM4, float);                 //! PM4 ZNA         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNAtdc, ZNATDC, float);                 //! TDC ZNA         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpmc, commonPMZNC, float);            //! PMC ZNC         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm1, ZNCPM1, float);                 //! PM1 ZNC         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm2, ZNCPM2, float);                 //! PM2 ZNC         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm3, ZNCPM3, float);                 //! PM3 ZNC         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm4, ZNCPM4, float);                 //! PM4 ZNC         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCtdc, ZNCTDC, float);                 //! TDC ZNC         // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(Centrality, centrality, float);         //! Centrality
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);        //! Timestamp
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t); //! Selection Flags
} // namespace zdcextra

DECLARE_SOA_TABLE(ZdcExtras, "AOD", "ZDCEXTRA", o2::soa::Index<>,
                  zdcextra::ZNApmc,
                  zdcextra::ZNApm1,
                  zdcextra::ZNApm2,
                  zdcextra::ZNApm3,
                  zdcextra::ZNApm4,
                  zdcextra::ZNAtdc,
                  zdcextra::ZNCpmc,
                  zdcextra::ZNCpm1,
                  zdcextra::ZNCpm2,
                  zdcextra::ZNCpm3,
                  zdcextra::ZNCpm4,
                  zdcextra::ZNCtdc,
                  zdcextra::Centrality,
                  zdcextra::Timestamp,
                  zdcextra::SelectionBits);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_ZDCEXTRA_H_
