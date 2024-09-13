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
#ifndef ZDC_INTERCALIB_H_
#define ZDC_INTERCALIB_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace vars
{
DECLARE_SOA_COLUMN(pmcZNA, ZNAcommonPM, float);      //! PMC ZNA
DECLARE_SOA_COLUMN(pm1ZNA, ZNAPM1, float);           //! PM1 ZNA
DECLARE_SOA_COLUMN(pm2ZNA, ZNAPM2, float);           //! PM2 ZNA
DECLARE_SOA_COLUMN(pm3ZNA, ZNAPM3, float);           //! PM3 ZNA
DECLARE_SOA_COLUMN(pm4ZNA, ZNAPM4, float);           //! PM4 ZNA
DECLARE_SOA_COLUMN(pmcZNC, ZNCcommonPM, float);      //! PMC ZNC
DECLARE_SOA_COLUMN(pm1ZNC, ZNCPM1, float);           //! PM1 ZNC
DECLARE_SOA_COLUMN(pm2ZNC, ZNCPM2, float);           //! PM2 ZNC
DECLARE_SOA_COLUMN(pm3ZNC, ZNCPM3, float);           //! PM3 ZNC
DECLARE_SOA_COLUMN(pm4ZNC, ZNCPM4, float);           //! PM4 ZNC

} // namespace vars

DECLARE_SOA_TABLE(ZDCInterCalib, "AOD", "ZDCIC", o2::soa::Index<>,
                  vars::pmcZNA,
                  vars::pm1ZNA,
                  vars::pm2ZNA,
                  vars::pm3ZNA,
                  vars::pm4ZNA,
                  vars::pmcZNC,
                  vars::pm1ZNC,
                  vars::pm2ZNC,
                  vars::pm3ZNC,
                  vars::pm4ZNC);
} // namespace o2::aod

#endif 
