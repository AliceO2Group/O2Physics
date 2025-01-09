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

#ifndef COMMON_DATAMODEL_ZDCINTERCALIB_H_
#define COMMON_DATAMODEL_ZDCINTERCALIB_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace znoutput
{
DECLARE_SOA_COLUMN(ZNApmc, commonPMZNA, float); //! PMC ZNA // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm1, ZNAPM1, float);      //! PM1 ZNA // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm2, ZNAPM2, float);      //! PM2 ZNA // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm3, ZNAPM3, float);      //! PM3 ZNA // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNApm4, ZNAPM4, float);      //! PM4 ZNA // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNAtdc, ZNATDC, float);      //! TDC ZNA // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpmc, commonPMZNC, float); //! PMC ZNC // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm1, ZNCPM1, float);      //! PM1 ZNC // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm2, ZNCPM2, float);      //! PM2 ZNC // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm3, ZNCPM3, float);      //! PM3 ZNC // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCpm4, ZNCPM4, float);      //! PM4 ZNC // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZNCtdc, ZNCTDC, float);      //! TDC ZNC // o2-linter: disable=name/o2-column

} // namespace znoutput

DECLARE_SOA_TABLE(ZDCInterCalib, "AOD", "ZDCIC", o2::soa::Index<>,
                  znoutput::ZNApmc,
                  znoutput::ZNApm1,
                  znoutput::ZNApm2,
                  znoutput::ZNApm3,
                  znoutput::ZNApm4,
                  znoutput::ZNAtdc,
                  znoutput::ZNCpmc,
                  znoutput::ZNCpm1,
                  znoutput::ZNCpm2,
                  znoutput::ZNCpm3,
                  znoutput::ZNCpm4,
                  znoutput::ZNCtdc);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_ZDCINTERCALIB_H_
