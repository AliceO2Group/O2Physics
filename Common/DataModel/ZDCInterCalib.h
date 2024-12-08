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
///
/// \brief Data model for ZN tower intercalibration
/// \author chiara.oppedisano@cern.ch
#ifndef COMMON_DATAMODEL_ZDCINTERCALIB_H_
#define COMMON_DATAMODEL_ZDCINTERCALIB_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace znoutput
{
DECLARE_SOA_COLUMN(ZNAtdc, znatdc, float); //! ZNA tdc value
DECLARE_SOA_COLUMN(ZNCtdc, znctdc, float); //! ZNC tdc value
DECLARE_SOA_COLUMN(ZNApmc, znapmc, float); //! PMC ZNA
DECLARE_SOA_COLUMN(ZNApm1, znapm1, float); //! PM1 ZNA
DECLARE_SOA_COLUMN(ZNApm2, znapm2, float); //! PM2 ZNA
DECLARE_SOA_COLUMN(ZNApm3, znapm3, float); //! PM3 ZNA
DECLARE_SOA_COLUMN(ZNApm4, znapm4, float); //! PM4 ZNA
DECLARE_SOA_COLUMN(ZNCpmc, zncpmc, float); //! PMC ZNC
DECLARE_SOA_COLUMN(ZNCpm1, zncpm1, float); //! PM1 ZNC
DECLARE_SOA_COLUMN(ZNCpm2, zncpm2, float); //! PM2 ZNC
DECLARE_SOA_COLUMN(ZNCpm3, zncpm3, float); //! PM3 ZNC
DECLARE_SOA_COLUMN(ZNCpm4, zncpm4, float); //! PM4 ZNC

} // namespace znoutput

DECLARE_SOA_TABLE(ZDCInterCalib, "AOD", "ZDCIC", o2::soa::Index<>,
                  znoutput::ZNAtdc,
                  znoutput::ZNCtdc,
                  znoutput::ZNApmc,
                  znoutput::ZNApm1,
                  znoutput::ZNApm2,
                  znoutput::ZNApm3,
                  znoutput::ZNApm4,
                  znoutput::ZNCpmc,
                  znoutput::ZNCpm1,
                  znoutput::ZNCpm2,
                  znoutput::ZNCpm3,
                  znoutput::ZNCpm4);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_ZDCINTERCALIB_H_
