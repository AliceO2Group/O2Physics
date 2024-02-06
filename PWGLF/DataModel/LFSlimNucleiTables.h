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
#include "Common/DataModel/Centrality.h"

#ifndef PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_

namespace o2::aod
{
namespace NucleiTableNS
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);
DECLARE_SOA_COLUMN(Beta, beta, float);
DECLARE_SOA_COLUMN(Zvertex, zVertex, float);
DECLARE_SOA_COLUMN(DCAxy, dcaxy, float);
DECLARE_SOA_COLUMN(DCAz, dcaz, float);
DECLARE_SOA_COLUMN(TPCsignal, tpcSignal, float);
DECLARE_SOA_COLUMN(ITSchi2, itsChi2, float);
DECLARE_SOA_COLUMN(TPCchi2, tpcChi2, float);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(TPCfindableCls, tpcFindableCls, uint8_t);
DECLARE_SOA_COLUMN(TPCcrossedRows, tpcCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(ITSclsMap, itsClsMap, uint8_t);
DECLARE_SOA_COLUMN(TPCnCls, tpcNCls, uint8_t);
DECLARE_SOA_COLUMN(ITSclusterSizes, itsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(gPt, genPt, float);
DECLARE_SOA_COLUMN(gEta, genEta, float);
DECLARE_SOA_COLUMN(gPhi, genPhi, float);
DECLARE_SOA_COLUMN(PDGcode, pdgCode, int);
DECLARE_SOA_COLUMN(SurvivedEventSelection, survivedEventSelection, bool);

} // namespace NucleiTableNS
namespace NucleiFlowTableNS
{
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(XQvecFV0A, xQvecFV0A, float);
DECLARE_SOA_COLUMN(YQvecFV0A, yQvecFV0A, float);
DECLARE_SOA_COLUMN(AmplQvecFV0A, amplQvecFV0A, float);
DECLARE_SOA_COLUMN(XQvecFT0M, xQvecFT0M, float);
DECLARE_SOA_COLUMN(YQvecFT0M, yQvecFT0M, float);
DECLARE_SOA_COLUMN(AmplQvecFT0M, amplQvecFT0M, float);
DECLARE_SOA_COLUMN(XQvecFT0A, xQvecFT0A, float);
DECLARE_SOA_COLUMN(YQvecFT0A, yQvecFT0A, float);
DECLARE_SOA_COLUMN(AmplQvecFT0A, amplQvecFT0A, float);
DECLARE_SOA_COLUMN(XQvecFT0C, xQvecFT0C, float);
DECLARE_SOA_COLUMN(YQvecFT0C, yQvecFT0C, float);
DECLARE_SOA_COLUMN(AmplQvecFT0C, amplQvecFT0C, float);
DECLARE_SOA_COLUMN(XQvecTPCpos, xQvecTPCpos, float);
DECLARE_SOA_COLUMN(YQvecTPCpos, yQvecTPCpos, float);
DECLARE_SOA_COLUMN(AmplQvecTPCpos, amplQvecTPCpos, float);
DECLARE_SOA_COLUMN(XQvecTPCneg, xQvecTPCneg, float);
DECLARE_SOA_COLUMN(YQvecTPCneg, yQvecTPCneg, float);
DECLARE_SOA_COLUMN(AmplQvecTPCneg, amplQvecTPCneg, float);
} // namespace NucleiFlowTableNS

DECLARE_SOA_TABLE(NucleiFlowColls, "AOD", "NUCLEIFLOWCOLL",
                  o2::soa::Index<>,
                  NucleiFlowTableNS::CentFV0A,
                  NucleiFlowTableNS::CentFT0M,
                  NucleiFlowTableNS::CentFT0A,
                  NucleiFlowTableNS::CentFT0C,
                  NucleiFlowTableNS::XQvecFV0A,
                  NucleiFlowTableNS::YQvecFV0A,
                  NucleiFlowTableNS::AmplQvecFV0A,
                  NucleiFlowTableNS::XQvecFT0M,
                  NucleiFlowTableNS::YQvecFT0M,
                  NucleiFlowTableNS::AmplQvecFT0M,
                  NucleiFlowTableNS::XQvecFT0A,
                  NucleiFlowTableNS::YQvecFT0A,
                  NucleiFlowTableNS::AmplQvecFT0A,
                  NucleiFlowTableNS::XQvecFT0C,
                  NucleiFlowTableNS::YQvecFT0C,
                  NucleiFlowTableNS::AmplQvecFT0C,
                  NucleiFlowTableNS::XQvecTPCpos,
                  NucleiFlowTableNS::YQvecTPCpos,
                  NucleiFlowTableNS::AmplQvecTPCpos,
                  NucleiFlowTableNS::XQvecTPCneg,
                  NucleiFlowTableNS::YQvecTPCneg,
                  NucleiFlowTableNS::AmplQvecTPCneg)

using NucleiFlowColl = NucleiFlowColls::iterator;

namespace NucleiTableNS
{
DECLARE_SOA_INDEX_COLUMN(NucleiFlowColl, nucleiFlowColl);
}

DECLARE_SOA_TABLE(NucleiTable, "AOD", "NUCLEITABLE",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::Phi,
                  NucleiTableNS::TPCInnerParam,
                  NucleiTableNS::Beta,
                  NucleiTableNS::Zvertex,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::TPCsignal,
                  NucleiTableNS::ITSchi2,
                  NucleiTableNS::TPCchi2,
                  NucleiTableNS::Flags,
                  NucleiTableNS::TPCfindableCls,
                  NucleiTableNS::TPCcrossedRows,
                  NucleiTableNS::ITSclsMap,
                  NucleiTableNS::TPCnCls,
                  NucleiTableNS::ITSclusterSizes,
                  NucleiTableNS::NucleiFlowCollId)

DECLARE_SOA_TABLE(NucleiTableMC, "AOD", "NUCLEITABLEMC",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::Phi,
                  NucleiTableNS::TPCInnerParam,
                  NucleiTableNS::Beta,
                  NucleiTableNS::Zvertex,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::TPCsignal,
                  NucleiTableNS::ITSchi2,
                  NucleiTableNS::TPCchi2,
                  NucleiTableNS::Flags,
                  NucleiTableNS::TPCfindableCls,
                  NucleiTableNS::TPCcrossedRows,
                  NucleiTableNS::ITSclsMap,
                  NucleiTableNS::TPCnCls,
                  NucleiTableNS::ITSclusterSizes,
                  NucleiTableNS::gPt,
                  NucleiTableNS::gEta,
                  NucleiTableNS::gPhi,
                  NucleiTableNS::PDGcode,
                  NucleiTableNS::SurvivedEventSelection)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_
