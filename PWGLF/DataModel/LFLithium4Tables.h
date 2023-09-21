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
/// \file LFLithium4Tables.h
/// \brief Slim tables for Lithium4
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFLITHIUM4TABLES_H_
#define PWGLF_DATAMODEL_LFLITHIUM4TABLES_H_

namespace o2::aod
{
namespace Lithium4TablesNS
{

DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);

DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);
DECLARE_SOA_COLUMN(PhiHe3, phiHe3, float);
DECLARE_SOA_COLUMN(PtPr, ptPr, float);
DECLARE_SOA_COLUMN(EtaPr, etaPr, float);
DECLARE_SOA_COLUMN(PhiPr, phiPr, float);

DECLARE_SOA_COLUMN(He3DCAxy, he3DcaXY, float);
DECLARE_SOA_COLUMN(He3DCAz, he3DcaZ, float);
DECLARE_SOA_COLUMN(ProtonDCAxy, protonDcaXY, float);
DECLARE_SOA_COLUMN(ProtonDCAz, protonDcaZ, float);

DECLARE_SOA_COLUMN(He3SignalTPC, he3SignalTPC, float);
DECLARE_SOA_COLUMN(He3InnerParamTPC, he3InnerParamTPC, float);
DECLARE_SOA_COLUMN(He3NClsTPC, he3NClsTPC, uint8_t);
DECLARE_SOA_COLUMN(He3NSigmaTPC, he3NSigmaTPC, float);
DECLARE_SOA_COLUMN(PrNSigmaTPC, prNSigmaTOF, float);
DECLARE_SOA_COLUMN(He3MassTOF, he3MassTOF, float);
DECLARE_SOA_COLUMN(PrMassTOF, prMassTOF, float);

DECLARE_SOA_COLUMN(IsBkgLS, isBkgLS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(SignedPtMC, signedPtMC, float);
DECLARE_SOA_COLUMN(MassMC, massMC, float);

} // namespace Lithium4TablesNS

DECLARE_SOA_TABLE(Lithium4Table, "AOD", "LITHIUM4TABLE",
                  Lithium4TablesNS::IsMatter,
                  Lithium4TablesNS::PtHe3,
                  Lithium4TablesNS::EtaHe3,
                  Lithium4TablesNS::PhiHe3,
                  Lithium4TablesNS::PtPr,
                  Lithium4TablesNS::EtaPr,
                  Lithium4TablesNS::PhiPr,
                  Lithium4TablesNS::He3DCAxy,
                  Lithium4TablesNS::He3DCAz,
                  Lithium4TablesNS::ProtonDCAxy,
                  Lithium4TablesNS::ProtonDCAz,
                  Lithium4TablesNS::He3SignalTPC,
                  Lithium4TablesNS::He3InnerParamTPC,
                  Lithium4TablesNS::He3NClsTPC,
                  Lithium4TablesNS::He3NSigmaTPC,
                  Lithium4TablesNS::PrNSigmaTPC,
                  Lithium4TablesNS::He3MassTOF,
                  Lithium4TablesNS::PrMassTOF,
                  Lithium4TablesNS::IsBkgLS,
                  Lithium4TablesNS::IsBkgEM)
DECLARE_SOA_TABLE(Lithium4TableMC, "AOD", "LITHIUM4TABLEMC",
                  Lithium4TablesNS::IsMatter,
                  Lithium4TablesNS::PtHe3,
                  Lithium4TablesNS::EtaHe3,
                  Lithium4TablesNS::PhiHe3,
                  Lithium4TablesNS::PtPr,
                  Lithium4TablesNS::EtaPr,
                  Lithium4TablesNS::PhiPr,
                  Lithium4TablesNS::He3DCAxy,
                  Lithium4TablesNS::He3DCAz,
                  Lithium4TablesNS::ProtonDCAxy,
                  Lithium4TablesNS::ProtonDCAz,
                  Lithium4TablesNS::He3SignalTPC,
                  Lithium4TablesNS::He3InnerParamTPC,
                  Lithium4TablesNS::He3NClsTPC,
                  Lithium4TablesNS::He3NSigmaTPC,
                  Lithium4TablesNS::PrNSigmaTPC,
                  Lithium4TablesNS::He3MassTOF,
                  Lithium4TablesNS::PrMassTOF,
                  Lithium4TablesNS::IsBkgLS,
                  Lithium4TablesNS::IsBkgEM,
                  Lithium4TablesNS::SignedPtMC,
                  Lithium4TablesNS::MassMC)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFLITHIUM4TABLES_H_
