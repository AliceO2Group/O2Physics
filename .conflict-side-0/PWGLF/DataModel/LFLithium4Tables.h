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

DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);
DECLARE_SOA_COLUMN(PhiHe3, phiHe3, float);
DECLARE_SOA_COLUMN(PtPr, ptPr, float);
DECLARE_SOA_COLUMN(EtaPr, etaPr, float);
DECLARE_SOA_COLUMN(PhiPr, phiPr, float);

DECLARE_SOA_COLUMN(DCAxyHe3, dcaxyHe3, float);
DECLARE_SOA_COLUMN(DCAzHe3, dcazHe3, float);
DECLARE_SOA_COLUMN(DCAxyPr, dcaxyPr, float);
DECLARE_SOA_COLUMN(DCAzPr, dcazPr, float);

DECLARE_SOA_COLUMN(SignalTPCHe3, signalTPCHe3, float);
DECLARE_SOA_COLUMN(InnerParamTPCHe3, innerParamTPCHe3, float);
DECLARE_SOA_COLUMN(SignalTPCPr, signalTPCPr, float);
DECLARE_SOA_COLUMN(InnerParamTPCPr, innerParamTPCPr, float);
DECLARE_SOA_COLUMN(NClsTPCHe3, nClsTPCHe3, uint8_t);
DECLARE_SOA_COLUMN(NSigmaTPCHe3, nSigmaTPCHe3, float);
DECLARE_SOA_COLUMN(NSigmaTPCPr, nSigmaTOFPr, float);
DECLARE_SOA_COLUMN(Chi2TPCHe3, chi2TPCHe3, float);
DECLARE_SOA_COLUMN(Chi2TPCPr, chi2TPCPr, float);
DECLARE_SOA_COLUMN(MassTOFHe3, massTOFHe3, float);
DECLARE_SOA_COLUMN(MassTOFPr, massTOFPr, float);
DECLARE_SOA_COLUMN(PIDtrkHe3, pidTrkHe3, uint32_t);
DECLARE_SOA_COLUMN(PIDtrkPr, pidTrkPr, uint32_t);

DECLARE_SOA_COLUMN(ItsClusterSizeHe3, itsClusterSizeHe3, uint32_t);
DECLARE_SOA_COLUMN(ItsClusterSizePr, itsClusterSizePr, uint32_t);

DECLARE_SOA_COLUMN(SharedClustersHe3, sharedClustersHe3, uint8_t);
DECLARE_SOA_COLUMN(SharedClustersPr, sharedClustersPr, uint8_t);

DECLARE_SOA_COLUMN(IsBkgLS, isBkgLS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(PtMCHe3, ptMCHe3, float);
DECLARE_SOA_COLUMN(EtaMCHe3, etaMCHe3, float);
DECLARE_SOA_COLUMN(PhiMCHe3, phiMCHe3, float);
DECLARE_SOA_COLUMN(PtMCPr, ptMCPr, float);
DECLARE_SOA_COLUMN(EtaMCPr, etaMCPr, float);
DECLARE_SOA_COLUMN(PhiMCPr, phiMCPr, float);
DECLARE_SOA_COLUMN(SignedPtMC, signedPtMC, float);
DECLARE_SOA_COLUMN(MassMC, massMC, float);

DECLARE_SOA_COLUMN(Multiplicity, multiplicity, uint16_t);
DECLARE_SOA_COLUMN(CentralityFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);

} // namespace Lithium4TablesNS

DECLARE_SOA_TABLE(Lithium4Table, "AOD", "LITHIUM4TABLE",
                  Lithium4TablesNS::PtHe3,
                  Lithium4TablesNS::EtaHe3,
                  Lithium4TablesNS::PhiHe3,
                  Lithium4TablesNS::PtPr,
                  Lithium4TablesNS::EtaPr,
                  Lithium4TablesNS::PhiPr,
                  Lithium4TablesNS::DCAxyHe3,
                  Lithium4TablesNS::DCAzHe3,
                  Lithium4TablesNS::DCAxyPr,
                  Lithium4TablesNS::DCAzPr,
                  Lithium4TablesNS::SignalTPCHe3,
                  Lithium4TablesNS::InnerParamTPCHe3,
                  Lithium4TablesNS::SignalTPCPr,
                  Lithium4TablesNS::InnerParamTPCPr,
                  Lithium4TablesNS::NClsTPCHe3,
                  Lithium4TablesNS::NSigmaTPCHe3,
                  Lithium4TablesNS::NSigmaTPCPr,
                  Lithium4TablesNS::Chi2TPCHe3,
                  Lithium4TablesNS::Chi2TPCPr,
                  Lithium4TablesNS::MassTOFHe3,
                  Lithium4TablesNS::MassTOFPr,
                  Lithium4TablesNS::PIDtrkHe3,
                  Lithium4TablesNS::PIDtrkPr,
                  Lithium4TablesNS::ItsClusterSizeHe3,
                  Lithium4TablesNS::ItsClusterSizePr,
                  Lithium4TablesNS::SharedClustersHe3,
                  Lithium4TablesNS::SharedClustersPr,
                  Lithium4TablesNS::IsBkgLS,
                  Lithium4TablesNS::IsBkgEM)
DECLARE_SOA_TABLE(Lithium4TableMC, "AOD", "LITHIUM4TABLEMC",
                  Lithium4TablesNS::PtMCHe3,
                  Lithium4TablesNS::EtaMCHe3,
                  Lithium4TablesNS::PhiMCHe3,
                  Lithium4TablesNS::PtMCPr,
                  Lithium4TablesNS::EtaMCPr,
                  Lithium4TablesNS::PhiMCPr,
                  Lithium4TablesNS::SignedPtMC,
                  Lithium4TablesNS::MassMC)
DECLARE_SOA_TABLE(Lithium4Mult, "AOD", "LITHIUM4MULT",
                  Lithium4TablesNS::Multiplicity,
                  Lithium4TablesNS::CentralityFT0C,
                  Lithium4TablesNS::MultiplicityFT0C)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFLITHIUM4TABLES_H_
