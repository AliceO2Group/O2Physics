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
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float); // centrality with FT0A estimator
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float); // centrality with FT0A estimator
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float); // centrality with FT0C estimator
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float); // centrality with FT0M estimator
DECLARE_SOA_COLUMN(PsiFT0A, psiFT0A, float);   // Psi with FT0A estimator
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float); // Multiplicity with FT0A estimator
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);   // Psi with FT0C estimator
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float); // Multiplicity with FT0C estimator
DECLARE_SOA_COLUMN(PsiTPC, psiTPC, float);     // Psi with TPC estimator
DECLARE_SOA_COLUMN(PsiTPCl, psiTPCl, float);   // Psi with TPC estimator (left)
DECLARE_SOA_COLUMN(PsiTPCr, psiTPCr, float);   // Psi with TPC estimator (right)
DECLARE_SOA_COLUMN(MultTPC, multTPC, int);     // Multiplicity with TPC estimator
} // namespace NucleiFlowTableNS

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
                  NucleiTableNS::ITSclusterSizes);

DECLARE_SOA_TABLE(NucleiTableFlow, "AOD", "NUCLEITABLEFLOW",
                  NucleiFlowTableNS::CentFV0A,
                  NucleiFlowTableNS::CentFT0M,
                  NucleiFlowTableNS::CentFT0A,
                  NucleiFlowTableNS::CentFT0C,
                  NucleiFlowTableNS::PsiFT0A,
                  NucleiFlowTableNS::MultFT0A,
                  NucleiFlowTableNS::PsiFT0C,
                  NucleiFlowTableNS::MultFT0C,
                  NucleiFlowTableNS::PsiTPC,
                  NucleiFlowTableNS::PsiTPCl,
                  NucleiFlowTableNS::PsiTPCr,
                  NucleiFlowTableNS::MultTPC);

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
