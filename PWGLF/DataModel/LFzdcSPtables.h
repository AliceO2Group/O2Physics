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

#ifndef PWGLF_DATAMODEL_ZDCSPTABLES_H_
#define PWGLF_DATAMODEL_ZDCSPTABLES_H_

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Centrality.h"

namespace o2::aod
{
namespace ZdcSPTableNS
{
DECLARE_SOA_COLUMN(EnergyZNA0, energyZNA0, float);
DECLARE_SOA_COLUMN(EnergyZNA1, energyZNA1, float);
DECLARE_SOA_COLUMN(EnergyZNA2, energyZNA2, float);
DECLARE_SOA_COLUMN(EnergyZNA3, energyZNA3, float);
DECLARE_SOA_COLUMN(EnergyZNA4, energyZNA4, float);
DECLARE_SOA_COLUMN(EnergyZNC0, energyZNC0, float);
DECLARE_SOA_COLUMN(EnergyZNC1, energyZNC1, float);
DECLARE_SOA_COLUMN(EnergyZNC2, energyZNC2, float);
DECLARE_SOA_COLUMN(EnergyZNC3, energyZNC3, float);
DECLARE_SOA_COLUMN(EnergyZNC4, energyZNC4, float);

} // namespace ZdcSPTableNS
DECLARE_SOA_TABLE(ZdcSPTable, "AOD", "ZDCSPTABLE",
                  bc::GlobalBC,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  cent::CentFT0C,
                  ZdcSPTableNS::EnergyZNA0,
                  ZdcSPTableNS::EnergyZNA1,
                  ZdcSPTableNS::EnergyZNA2,
                  ZdcSPTableNS::EnergyZNA3,
                  ZdcSPTableNS::EnergyZNA4,
                  ZdcSPTableNS::EnergyZNC0,
                  ZdcSPTableNS::EnergyZNC1,
                  ZdcSPTableNS::EnergyZNC2,
                  ZdcSPTableNS::EnergyZNC3,
                  ZdcSPTableNS::EnergyZNC4);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_ZDCSPTABLES_H_
