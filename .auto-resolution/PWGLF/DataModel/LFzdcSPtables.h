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

#ifndef PWGLF_DATAMODEL_LFZDCSPTABLES_H_
#define PWGLF_DATAMODEL_LFZDCSPTABLES_H_

#include "Common/DataModel/Centrality.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace ZdcSPTableNS
{
DECLARE_SOA_COLUMN(TimeSinceSOR, timeSinceSOR, uint64_t);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, float);
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
DECLARE_SOA_COLUMN(AmpZNA0, ampZNA0, float);
DECLARE_SOA_COLUMN(AmpZNA1, ampZNA1, float);
DECLARE_SOA_COLUMN(AmpZNA2, ampZNA2, float);
DECLARE_SOA_COLUMN(AmpZNA3, ampZNA3, float);
DECLARE_SOA_COLUMN(AmpZNA4, ampZNA4, float);
DECLARE_SOA_COLUMN(AmpZNC0, ampZNC0, float);
DECLARE_SOA_COLUMN(AmpZNC1, ampZNC1, float);
DECLARE_SOA_COLUMN(AmpZNC2, ampZNC2, float);
DECLARE_SOA_COLUMN(AmpZNC3, ampZNC3, float);
DECLARE_SOA_COLUMN(AmpZNC4, ampZNC4, float);

} // namespace ZdcSPTableNS
DECLARE_SOA_TABLE(ZdcSPTable, "AOD", "ZDCSPTABLE",
                  ZdcSPTableNS::TimeSinceSOR,
                  bc::RunNumber,
                  ZdcSPTableNS::HadronicRate,
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

#endif // PWGLF_DATAMODEL_LFZDCSPTABLES_H_
