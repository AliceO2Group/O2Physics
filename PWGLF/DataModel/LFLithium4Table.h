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
/// \file LFLithium4Table.h
/// \brief Slim table for Lithium4
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFLITHIUM4TABLE_H_
#define PWGLF_DATAMODEL_LFLITHIUM4TABLE_H_

namespace o2::aod
{
namespace Lithium4TableNS
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rapidity, rapidity, float);
DECLARE_SOA_COLUMN(Mass, mass, float);

DECLARE_SOA_COLUMN(He3DCAxy, he3DcaXY, float);
DECLARE_SOA_COLUMN(He3DCAz, he3DcaZ, float);
DECLARE_SOA_COLUMN(ProtonDCAxy, protonDcaXY, float);
DECLARE_SOA_COLUMN(ProtonDCAz, protonDcaZ, float);

DECLARE_SOA_COLUMN(He3SignalTPC, he3SignalTPC, float);
DECLARE_SOA_COLUMN(He3InnerParamTPC, he3InnerParamTPC, float);
DECLARE_SOA_COLUMN(He3NClsTPC, he3NClsTPC, uint8_t);
DECLARE_SOA_COLUMN(He3NSigmaTPC, he3NSigmaTPC, float);

DECLARE_SOA_COLUMN(IsBkgLS, isBkgLS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(PtMC, ptMC, float);
DECLARE_SOA_COLUMN(MassMC, massMC, float);

} // namespace Lithium4TableNS

DECLARE_SOA_TABLE(Lithium4Table, "AOD", "LITHIUM4TABLE",
                  Lithium4TableNS::Pt,
                  Lithium4TableNS::Rapidity,
                  Lithium4TableNS::Mass,
                  Lithium4TableNS::He3DCAxy,
                  Lithium4TableNS::He3DCAz,
                  Lithium4TableNS::ProtonDCAxy,
                  Lithium4TableNS::ProtonDCAz,
                  Lithium4TableNS::He3SignalTPC,
                  Lithium4TableNS::He3InnerParamTPC,
                  Lithium4TableNS::He3NClsTPC,
                  Lithium4TableNS::He3NSigmaTPC,
                  Lithium4TableNS::IsBkgLS,
                  Lithium4TableNS::IsBkgEM)
DECLARE_SOA_TABLE(Lithium4TableMC, "AOD", "LITHIUM4TABLEMC",
                  Lithium4TableNS::Pt,
                  Lithium4TableNS::Rapidity,
                  Lithium4TableNS::Mass,
                  Lithium4TableNS::He3DCAxy,
                  Lithium4TableNS::He3DCAz,
                  Lithium4TableNS::ProtonDCAxy,
                  Lithium4TableNS::ProtonDCAz,
                  Lithium4TableNS::He3SignalTPC,
                  Lithium4TableNS::He3InnerParamTPC,
                  Lithium4TableNS::He3NClsTPC,
                  Lithium4TableNS::He3NSigmaTPC,
                  Lithium4TableNS::IsBkgLS,
                  Lithium4TableNS::IsBkgEM,
                  Lithium4TableNS::PtMC,
                  Lithium4TableNS::MassMC)
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFLITHIUM4TABLE_H_
