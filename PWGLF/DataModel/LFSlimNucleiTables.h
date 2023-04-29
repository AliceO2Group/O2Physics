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
/// \file LFSlimNucleiTables.h
/// \brief Slim nuclei tables
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_

namespace o2::aod
{
namespace NucleiTableNS
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(ITSclsMap, itsClsMap, uint8_t);
DECLARE_SOA_COLUMN(TPCnCls, tpcNCls, uint8_t);
DECLARE_SOA_COLUMN(DCAxy, dcaxy, int8_t);
DECLARE_SOA_COLUMN(DCAz, dcaz, int8_t);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(TPCnsigma, tpcnsigma, uint8_t);
DECLARE_SOA_COLUMN(TOFmass, tofmass, uint8_t);
DECLARE_SOA_COLUMN(gPt, genPt, float);
DECLARE_SOA_COLUMN(gEta, genEta, float);
DECLARE_SOA_COLUMN(PDGcode, pdgCode, int);

} // namespace NucleiTableNS
DECLARE_SOA_TABLE(NucleiTable, "AOD", "NUCLEITABLE",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::ITSclsMap,
                  NucleiTableNS::TPCnCls,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::Flags,
                  NucleiTableNS::TPCnsigma,
                  NucleiTableNS::TOFmass)

DECLARE_SOA_TABLE(NucleiTableMC, "AOD", "NUCLEITABLEMC",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::ITSclsMap,
                  NucleiTableNS::TPCnCls,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::Flags,
                  NucleiTableNS::TPCnsigma,
                  NucleiTableNS::TOFmass,
                  NucleiTableNS::gPt,
                  NucleiTableNS::gEta,
                  NucleiTableNS::PDGcode)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_
