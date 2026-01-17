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
/// \file LFhe3HadronTables.h
/// \brief Slim tables for he3Hadron
///

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFHE3HADRONTABLES_H_
#define PWGLF_DATAMODEL_LFHE3HADRONTABLES_H_

namespace o2::aod
{
namespace he3HadronTablesNS
{

DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);
DECLARE_SOA_COLUMN(PhiHe3, phiHe3, float);
DECLARE_SOA_COLUMN(PtHad, ptHad, float);
DECLARE_SOA_COLUMN(EtaHad, etaHad, float);
DECLARE_SOA_COLUMN(PhiHad, phiHad, float);

DECLARE_SOA_COLUMN(DCAxyHe3, dcaxyHe3, float);
DECLARE_SOA_COLUMN(DCAzHe3, dcazHe3, float);
DECLARE_SOA_COLUMN(DCAxyHad, dcaxyHad, float);
DECLARE_SOA_COLUMN(DCAzHad, dcazHad, float);
DECLARE_SOA_COLUMN(DCApair, dcapair, float);

DECLARE_SOA_COLUMN(SignalTPCHe3, signalTPCHe3, float);
DECLARE_SOA_COLUMN(InnerParamTPCHe3, innerParamTPCHe3, float);
DECLARE_SOA_COLUMN(SignalTPCHad, signalTPCHad, float);
DECLARE_SOA_COLUMN(InnerParamTPCHad, innerParamTPCHad, float);
DECLARE_SOA_COLUMN(NClsTPCHe3, nClsTPCHe3, uint8_t);
DECLARE_SOA_COLUMN(NSigmaTPCHe3, nSigmaTPCHe3, float);
DECLARE_SOA_COLUMN(NSigmaTPCHadPi, nSigmaTPCHadPi, float);
DECLARE_SOA_COLUMN(NSigmaTPCHadKa, nSigmaTPCHadKa, float);
DECLARE_SOA_COLUMN(NSigmaTPCHadPr, nSigmaTPCHadPr, float);
DECLARE_SOA_COLUMN(NSigmaTOFHadPi, nSigmaTOFHadPi, float);
DECLARE_SOA_COLUMN(NSigmaTOFHadKa, nSigmaTOFHadKa, float);
DECLARE_SOA_COLUMN(NSigmaTOFHadPr, nSigmaTOFHadPr, float);
DECLARE_SOA_COLUMN(Chi2TPCHe3, chi2TPCHe3, float);
DECLARE_SOA_COLUMN(Chi2TPCHad, chi2TPCHad, float);
DECLARE_SOA_COLUMN(MassTOFHe3, massTOFHe3, float);
DECLARE_SOA_COLUMN(MassTOFHad, massTOFHad, float);
DECLARE_SOA_COLUMN(PIDtrkHe3, pidTrkHe3, uint32_t);
DECLARE_SOA_COLUMN(PIDtrkHad, pidTrkHad, uint32_t);
DECLARE_SOA_COLUMN(TrackIDHe3, trackIDHe3, int);
DECLARE_SOA_COLUMN(TrackIDHad, trackIDHad, int);

DECLARE_SOA_COLUMN(ItsClusterSizeHe3, itsClusterSizeHe3, uint32_t);
DECLARE_SOA_COLUMN(ItsClusterSizeHad, itsClusterSizeHad, uint32_t);

DECLARE_SOA_COLUMN(SharedClustersHe3, sharedClustersHe3, uint8_t);
DECLARE_SOA_COLUMN(SharedClustersHad, sharedClustersHad, uint8_t);

DECLARE_SOA_COLUMN(IsBkgUS, isBkgUS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(PtMCHe3, ptMCHe3, float);
DECLARE_SOA_COLUMN(EtaMCHe3, etaMCHe3, float);
DECLARE_SOA_COLUMN(PhiMCHe3, phiMCHe3, float);
DECLARE_SOA_COLUMN(PtMCHad, ptMCHad, float);
DECLARE_SOA_COLUMN(EtaMCHad, etaMCHad, float);
DECLARE_SOA_COLUMN(PhiMCHad, phiMCHad, float);
DECLARE_SOA_COLUMN(SignedPtMC, signedPtMC, float);
DECLARE_SOA_COLUMN(MassMC, massMC, float);

DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(ZVertex, zVertex, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, uint16_t);
DECLARE_SOA_COLUMN(CentralityFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);

/* Flags: 0 - both primary,
          1 - both from Li4,
          2 - both from hypertriton,
          3 - mixed pair (a primary and one from Li4/hypertriton/material/other decays or any other combination)
*/
DECLARE_SOA_COLUMN(Flags, flags, uint8_t);

} // namespace he3HadronTablesNS

DECLARE_SOA_TABLE(he3HadronTable, "AOD", "HE3HADTABLE",
                  he3HadronTablesNS::PtHe3,
                  he3HadronTablesNS::EtaHe3,
                  he3HadronTablesNS::PhiHe3,
                  he3HadronTablesNS::PtHad,
                  he3HadronTablesNS::EtaHad,
                  he3HadronTablesNS::PhiHad,
                  he3HadronTablesNS::DCAxyHe3,
                  he3HadronTablesNS::DCAzHe3,
                  he3HadronTablesNS::DCAxyHad,
                  he3HadronTablesNS::DCAzHad,
                  he3HadronTablesNS::DCApair,
                  he3HadronTablesNS::SignalTPCHe3,
                  he3HadronTablesNS::InnerParamTPCHe3,
                  he3HadronTablesNS::SignalTPCHad,
                  he3HadronTablesNS::InnerParamTPCHad,
                  he3HadronTablesNS::NClsTPCHe3,
                  he3HadronTablesNS::NSigmaTPCHe3,
                  he3HadronTablesNS::NSigmaTPCHadPi,
                  he3HadronTablesNS::NSigmaTPCHadKa,
                  he3HadronTablesNS::NSigmaTPCHadPr,
                  he3HadronTablesNS::NSigmaTOFHadPi,
                  he3HadronTablesNS::NSigmaTOFHadKa,
                  he3HadronTablesNS::NSigmaTOFHadPr,
                  he3HadronTablesNS::Chi2TPCHe3,
                  he3HadronTablesNS::Chi2TPCHad,
                  he3HadronTablesNS::MassTOFHe3,
                  he3HadronTablesNS::MassTOFHad,
                  he3HadronTablesNS::PIDtrkHe3,
                  he3HadronTablesNS::PIDtrkHad,
                  he3HadronTablesNS::ItsClusterSizeHe3,
                  he3HadronTablesNS::ItsClusterSizeHad,
                  he3HadronTablesNS::SharedClustersHe3,
                  he3HadronTablesNS::SharedClustersHad)
DECLARE_SOA_TABLE(he3HadronTableMC, "AOD", "HE3HADTABLEMC",
                  he3HadronTablesNS::PtMCHe3,
                  he3HadronTablesNS::EtaMCHe3,
                  he3HadronTablesNS::PhiMCHe3,
                  he3HadronTablesNS::PtMCHad,
                  he3HadronTablesNS::EtaMCHad,
                  he3HadronTablesNS::PhiMCHad,
                  he3HadronTablesNS::SignedPtMC,
                  he3HadronTablesNS::MassMC,
                  he3HadronTablesNS::Flags)
DECLARE_SOA_TABLE(he3HadronMult, "AOD", "HE3HADMULT",
                  he3HadronTablesNS::CollisionId,
                  he3HadronTablesNS::ZVertex,
                  he3HadronTablesNS::Multiplicity,
                  he3HadronTablesNS::CentralityFT0C,
                  he3HadronTablesNS::MultiplicityFT0C)
DECLARE_SOA_TABLE(he3HadronQa, "AOD", "HE3HADQA",
                  he3HadronTablesNS::TrackIDHe3,
                  he3HadronTablesNS::TrackIDHad)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFHE3HADRONTABLES_H_
