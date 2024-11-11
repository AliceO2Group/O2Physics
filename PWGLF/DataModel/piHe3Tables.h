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
/// \file PiHe3Tables.h
/// \brief Slim tables for PiHe3
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFPIHE3TABLES_H_
#define PWGLF_DATAMODEL_LFPIHE3TABLES_H_

namespace o2::aod
{
namespace PiHe3TablesNS
{

DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);
DECLARE_SOA_COLUMN(PhiHe3, phiHe3, float);
DECLARE_SOA_COLUMN(PtPion, ptPion, float);
DECLARE_SOA_COLUMN(EtaPion, etaPion, float);
DECLARE_SOA_COLUMN(PhiPion, phiPion, float);

DECLARE_SOA_COLUMN(DCAxyHe3, dcaxyHe3, float);
DECLARE_SOA_COLUMN(DCAzHe3, dcazHe3, float);
DECLARE_SOA_COLUMN(DCAxyPion, dcaxyPion, float);
DECLARE_SOA_COLUMN(DCAzPion, dcazPion, float);

DECLARE_SOA_COLUMN(SignalTPCHe3, signalTPCHe3, float);
DECLARE_SOA_COLUMN(InnerParamTPCHe3, innerParamTPCHe3, float);
DECLARE_SOA_COLUMN(SignalTPCPion, signalTPCPion, float);
DECLARE_SOA_COLUMN(InnerParamTPCPion, innerParamTPCPion, float);
DECLARE_SOA_COLUMN(NClsTPCHe3, nClsTPCHe3, uint8_t);
DECLARE_SOA_COLUMN(NSigmaTPCHe3, nSigmaTPCHe3, float);
DECLARE_SOA_COLUMN(NSigmaTPCPion, nSigmaTOFPion, float);
DECLARE_SOA_COLUMN(Chi2TPCHe3, chi2TPCHe3, float);
DECLARE_SOA_COLUMN(Chi2TPCPion, chi2TPCPion, float);
DECLARE_SOA_COLUMN(MassTOFHe3, massTOFHe3, float);
DECLARE_SOA_COLUMN(MassTOFPion, massTOFPion, float);
DECLARE_SOA_COLUMN(PIDtrkHe3, pidTrkHe3, uint32_t);
DECLARE_SOA_COLUMN(PIDtrkPion, pidTrkPion, uint32_t);

DECLARE_SOA_COLUMN(ItsClusterSizeHe3, itsClusterSizeHe3, uint32_t);
DECLARE_SOA_COLUMN(ItsClusterSizePion, itsClusterSizePion, uint32_t);

DECLARE_SOA_COLUMN(SharedClustersHe3, sharedClustersHe3, uint8_t);
DECLARE_SOA_COLUMN(SharedClustersPion, sharedClustersPion, uint8_t);

DECLARE_SOA_COLUMN(IsBkgLS, isBkgLS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(PtMCHe3, ptMCHe3, float);
DECLARE_SOA_COLUMN(EtaMCHe3, etaMCHe3, float);
DECLARE_SOA_COLUMN(PhiMCHe3, phiMCHe3, float);
DECLARE_SOA_COLUMN(PtMCPion, ptMCPion, float);
DECLARE_SOA_COLUMN(EtaMCPion, etaMCPion, float);
DECLARE_SOA_COLUMN(PhiMCPion, phiMCPion, float);
DECLARE_SOA_COLUMN(SignedPtMC, signedPtMC, float);
DECLARE_SOA_COLUMN(MassMC, massMC, float);

DECLARE_SOA_COLUMN(Multiplicity, multiplicity, uint16_t);
DECLARE_SOA_COLUMN(CentralityFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);

} // namespace PiHe3TablesNS

DECLARE_SOA_TABLE(PiHe3Table, "AOD", "PIHE3TABLE",
                  PiHe3TablesNS::PtHe3,
                  PiHe3TablesNS::EtaHe3,
                  PiHe3TablesNS::PhiHe3,
                  PiHe3TablesNS::PtPion,
                  PiHe3TablesNS::EtaPion,
                  PiHe3TablesNS::PhiPion,
                  PiHe3TablesNS::DCAxyHe3,
                  PiHe3TablesNS::DCAzHe3,
                  PiHe3TablesNS::DCAxyPion,
                  PiHe3TablesNS::DCAzPion,
                  PiHe3TablesNS::SignalTPCHe3,
                  PiHe3TablesNS::InnerParamTPCHe3,
                  PiHe3TablesNS::SignalTPCPion,
                  PiHe3TablesNS::InnerParamTPCPion,
                  PiHe3TablesNS::NClsTPCHe3,
                  PiHe3TablesNS::NSigmaTPCHe3,
                  PiHe3TablesNS::NSigmaTPCPion,
                  PiHe3TablesNS::Chi2TPCHe3,
                  PiHe3TablesNS::Chi2TPCPion,
                  PiHe3TablesNS::MassTOFHe3,
                  PiHe3TablesNS::MassTOFPion,
                  PiHe3TablesNS::PIDtrkHe3,
                  PiHe3TablesNS::PIDtrkPion,
                  PiHe3TablesNS::ItsClusterSizeHe3,
                  PiHe3TablesNS::ItsClusterSizePion,
                  PiHe3TablesNS::SharedClustersHe3,
                  PiHe3TablesNS::SharedClustersPion,
                  PiHe3TablesNS::IsBkgLS,
                  PiHe3TablesNS::IsBkgEM)
DECLARE_SOA_TABLE(PiHe3TableMC, "AOD", "PIHE3TABLEMC",
                  PiHe3TablesNS::PtMCHe3,
                  PiHe3TablesNS::EtaMCHe3,
                  PiHe3TablesNS::PhiMCHe3,
                  PiHe3TablesNS::PtMCPion,
                  PiHe3TablesNS::EtaMCPion,
                  PiHe3TablesNS::PhiMCPion,
                  PiHe3TablesNS::SignedPtMC,
                  PiHe3TablesNS::MassMC)
DECLARE_SOA_TABLE(PiHe3Mult, "AOD", "PIHE3MULT",
                  PiHe3TablesNS::Multiplicity,
                  PiHe3TablesNS::CentralityFT0C,
                  PiHe3TablesNS::MultiplicityFT0C)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFPIHE3TABLES_H_
