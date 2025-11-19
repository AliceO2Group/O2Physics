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

#include "Common/DataModel/Centrality.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

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
DECLARE_SOA_COLUMN(NContrib, nContrib, int);
DECLARE_SOA_COLUMN(DCAxy, dcaxy, float);
DECLARE_SOA_COLUMN(DCAz, dcaz, float);
DECLARE_SOA_COLUMN(TPCsignal, tpcSignal, float);
DECLARE_SOA_COLUMN(ITSchi2, itsChi2, float);
DECLARE_SOA_COLUMN(TPCchi2, tpcChi2, float);
DECLARE_SOA_COLUMN(TOFchi2, tofChi2, float);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(TPCfindableCls, tpcFindableCls, uint8_t);
DECLARE_SOA_COLUMN(TPCcrossedRows, tpcCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(ITSclsMap, itsClsMap, uint8_t);
DECLARE_SOA_COLUMN(TPCnCls, tpcNCls, uint8_t);
DECLARE_SOA_COLUMN(TPCnClsShared, tpcNClsShared, uint8_t);
DECLARE_SOA_COLUMN(ITSclusterSizes, itsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(SurvivedEventSelection, survivedEventSelection, bool);
DECLARE_SOA_COLUMN(gPt, genPt, float);
DECLARE_SOA_COLUMN(gEta, genEta, float);
DECLARE_SOA_COLUMN(gPhi, genPhi, float);
DECLARE_SOA_COLUMN(PDGcode, pdgCode, int);
DECLARE_SOA_COLUMN(MotherPDGcode, MotherpdgCode, int);
DECLARE_SOA_COLUMN(MotherDecRad, motherDecRad, float);
DECLARE_SOA_COLUMN(AbsoDecL, absoDecL, float);
DECLARE_SOA_COLUMN(McProcess, mcProcess, uint64_t);

} // namespace NucleiTableNS

namespace NucleiPairTableNS
{
DECLARE_SOA_COLUMN(Pt1, pt1, float);                              // first particle pt
DECLARE_SOA_COLUMN(Eta1, eta1, float);                            // first particle eta
DECLARE_SOA_COLUMN(Phi1, phi1, float);                            // first particle phi
DECLARE_SOA_COLUMN(TPCInnerParam1, tpcInnerParam1, float);        // first particle TPC inner param
DECLARE_SOA_COLUMN(TPCsignal1, tpcSignal1, float);                // first particle TPC signal
DECLARE_SOA_COLUMN(DCAxy1, dcaxy1, float);                        // first particle DCA xy
DECLARE_SOA_COLUMN(DCAz1, dcaz1, float);                          // first particle DCA z
DECLARE_SOA_COLUMN(ClusterSizesITS1, clusterSizesITS1, uint32_t); // first particle ITS cluster sizes
DECLARE_SOA_COLUMN(Flags1, flags1, uint16_t);                     // first particle flags
DECLARE_SOA_COLUMN(Pt2, pt2, float);                              // second particle pt
DECLARE_SOA_COLUMN(Eta2, eta2, float);                            // second particle eta
DECLARE_SOA_COLUMN(Phi2, phi2, float);                            // second particle phi
DECLARE_SOA_COLUMN(TPCInnerParam2, tpcInnerParam2, float);        // second particle TPC inner param
DECLARE_SOA_COLUMN(TPCsignal2, tpcSignal2, float);                // second particle TPC signal
DECLARE_SOA_COLUMN(DCAxy2, dcaxy2, float);                        // second particle DCA xy
DECLARE_SOA_COLUMN(DCAz2, dcaz2, float);                          // second particle DCA z
DECLARE_SOA_COLUMN(ClusterSizesITS2, clusterSizesITS2, uint32_t); // second particle ITS cluster sizes
DECLARE_SOA_COLUMN(Flags2, flags2, uint16_t);                     // second particle flags
} // namespace NucleiPairTableNS

namespace NucleiFlowTableNS
{
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float); // centrality with FT0A estimator
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float); // centrality with FT0A estimator
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float); // centrality with FT0C estimator
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float); // centrality with FT0M estimator
DECLARE_SOA_COLUMN(PsiFT0A, psiFT0A, float);   // Psi with FT0A estimator
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);   // Psi with FT0C estimator
DECLARE_SOA_COLUMN(PsiTPC, psiTPC, float);     // Psi with TPC estimator
DECLARE_SOA_COLUMN(PsiTPCl, psiTPCl, float);   // Psi with TPC estimator (left)
DECLARE_SOA_COLUMN(PsiTPCr, psiTPCr, float);   // Psi with TPC estimator (right)
DECLARE_SOA_COLUMN(QFT0A, qFT0A, float);       // Amplitude with FT0A estimator
DECLARE_SOA_COLUMN(QFT0C, qFT0C, float);       // Amplitude with FT0C estimator
DECLARE_SOA_COLUMN(QTPC, qTPC, float);         // Amplitude with TPC estimator
DECLARE_SOA_COLUMN(QTPCl, qTPCl, float);       // Amplitude with TPC estimator (left)
DECLARE_SOA_COLUMN(QTPCr, qTPCr, float);       // Amplitude with TPC estimator (right)
} // namespace NucleiFlowTableNS

DECLARE_SOA_TABLE(NucleiTable, "AOD", "NUCLEITABLE",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::Phi,
                  NucleiTableNS::TPCInnerParam,
                  NucleiTableNS::Beta,
                  NucleiTableNS::Zvertex,
                  NucleiTableNS::NContrib,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::TPCsignal,
                  NucleiTableNS::ITSchi2,
                  NucleiTableNS::TPCchi2,
                  NucleiTableNS::TOFchi2,
                  NucleiTableNS::Flags,
                  NucleiTableNS::TPCfindableCls,
                  NucleiTableNS::TPCcrossedRows,
                  NucleiTableNS::ITSclsMap,
                  NucleiTableNS::TPCnCls,
                  NucleiTableNS::TPCnClsShared,
                  NucleiTableNS::ITSclusterSizes);

DECLARE_SOA_TABLE(NucleiTableFlow, "AOD", "NUCLEITABLEFLOW",
                  NucleiFlowTableNS::CentFV0A,
                  NucleiFlowTableNS::CentFT0M,
                  NucleiFlowTableNS::CentFT0A,
                  NucleiFlowTableNS::CentFT0C,
                  NucleiFlowTableNS::PsiFT0A,
                  NucleiFlowTableNS::PsiFT0C,
                  NucleiFlowTableNS::PsiTPC,
                  NucleiFlowTableNS::PsiTPCl,
                  NucleiFlowTableNS::PsiTPCr,
                  NucleiFlowTableNS::QFT0A,
                  NucleiFlowTableNS::QFT0C,
                  NucleiFlowTableNS::QTPC,
                  NucleiFlowTableNS::QTPCl,
                  NucleiFlowTableNS::QTPCr);

DECLARE_SOA_TABLE(NucleiTableMC, "AOD", "NUCLEITABLEMC",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::Phi,
                  NucleiTableNS::TPCInnerParam,
                  NucleiTableNS::Beta,
                  NucleiTableNS::Zvertex,
                  NucleiTableNS::NContrib,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::TPCsignal,
                  NucleiTableNS::ITSchi2,
                  NucleiTableNS::TPCchi2,
                  NucleiTableNS::TOFchi2,
                  NucleiTableNS::Flags,
                  NucleiTableNS::TPCfindableCls,
                  NucleiTableNS::TPCcrossedRows,
                  NucleiTableNS::ITSclsMap,
                  NucleiTableNS::TPCnCls,
                  NucleiTableNS::TPCnClsShared,
                  NucleiTableNS::ITSclusterSizes,
                  NucleiTableNS::SurvivedEventSelection,
                  NucleiTableNS::gPt,
                  NucleiTableNS::gEta,
                  NucleiTableNS::gPhi,
                  NucleiTableNS::PDGcode,
                  NucleiTableNS::MotherPDGcode,
                  NucleiTableNS::MotherDecRad,
                  NucleiTableNS::AbsoDecL);

DECLARE_SOA_TABLE(NucleiPairTable, "AOD", "NUCLEIPAIRTABLE",
                  NucleiPairTableNS::Pt1,
                  NucleiPairTableNS::Eta1,
                  NucleiPairTableNS::Phi1,
                  NucleiPairTableNS::TPCInnerParam1,
                  NucleiPairTableNS::TPCsignal1,
                  NucleiPairTableNS::DCAxy1,
                  NucleiPairTableNS::DCAz1,
                  NucleiPairTableNS::ClusterSizesITS1,
                  NucleiPairTableNS::Flags1,
                  NucleiPairTableNS::Pt2,
                  NucleiPairTableNS::Eta2,
                  NucleiPairTableNS::Phi2,
                  NucleiPairTableNS::TPCInnerParam2,
                  NucleiPairTableNS::TPCsignal2,
                  NucleiPairTableNS::DCAxy2,
                  NucleiPairTableNS::DCAz2,
                  NucleiPairTableNS::ClusterSizesITS2,
                  NucleiPairTableNS::Flags2);

// Reduced table
DECLARE_SOA_TABLE(NucleiTableRed, "AOD", "NUCLEITABLERED",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::Phi,
                  NucleiTableNS::TPCInnerParam,
                  NucleiTableNS::ITSclusterSizes,
                  NucleiTableNS::TPCsignal,
                  NucleiTableNS::Beta,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::Flags,
                  NucleiTableNS::McProcess,
                  NucleiTableNS::PDGcode,
                  NucleiTableNS::MotherPDGcode);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_
