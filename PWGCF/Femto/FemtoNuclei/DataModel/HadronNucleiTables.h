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
/// \file HadronNucleiTables.h
/// \brief Slim tables for piNuclei
/// \author CMY
/// \date 2025-04-10

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_HADRONNUCLEITABLES_H_
#define PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_HADRONNUCLEITABLES_H_

namespace o2::aod
{
namespace hadron_nuclei_tables
{

DECLARE_SOA_COLUMN(PtNu, ptNu, float);
DECLARE_SOA_COLUMN(EtaNu, etaNu, float);
DECLARE_SOA_COLUMN(PhiNu, phiNu, float);
DECLARE_SOA_COLUMN(PtHyp, ptHyp, float);
DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);
DECLARE_SOA_COLUMN(EtaHyp, etaHyp, float);
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);
DECLARE_SOA_COLUMN(PhiHyp, phiHyp, float);
DECLARE_SOA_COLUMN(PtHad, ptHad, float);
DECLARE_SOA_COLUMN(EtaHad, etaHad, float);
DECLARE_SOA_COLUMN(PhiHad, phiHad, float);

DECLARE_SOA_COLUMN(DcaxyNu, dcaxyNu, float);
DECLARE_SOA_COLUMN(DcazNu, dcazNu, float);
DECLARE_SOA_COLUMN(DcaxyHad, dcaxyHad, float);
DECLARE_SOA_COLUMN(DcazHad, dcazHad, float);

DECLARE_SOA_COLUMN(SignalTPCNu, signalTPCNu, float);
DECLARE_SOA_COLUMN(InnerParamTPCNu, innerParamTPCNu, float);
DECLARE_SOA_COLUMN(SignalTPCHad, signalTPCHad, float);
DECLARE_SOA_COLUMN(InnerParamTPCHad, innerParamTPCHad, float);
DECLARE_SOA_COLUMN(NClsTPCNu, nClsTPCNu, uint8_t);
DECLARE_SOA_COLUMN(NSigmaTPCNu, nSigmaTPCNu, float);
DECLARE_SOA_COLUMN(NSigmaTPCHad, nSigmaTPCHad, float);
DECLARE_SOA_COLUMN(Chi2TPCNu, chi2TPCNu, float);
DECLARE_SOA_COLUMN(Chi2TPCHad, chi2TPCHad, float);
DECLARE_SOA_COLUMN(MassTOFNu, massTOFNu, float);
DECLARE_SOA_COLUMN(MassTOFHad, massTOFHad, float);
DECLARE_SOA_COLUMN(HaddTrkNu, pidTrkNu, uint32_t);
DECLARE_SOA_COLUMN(HaddTrkHad, pidTrkHad, uint32_t);
DECLARE_SOA_COLUMN(TrackIDHad, trackIDHad, int);
DECLARE_SOA_COLUMN(TrackIDNu, trackIDNu, int);

DECLARE_SOA_COLUMN(ItsClusterSizeNu, itsClusterSizeNu, uint32_t);
DECLARE_SOA_COLUMN(ItsClusterSizeHad, itsClusterSizeHad, uint32_t);

DECLARE_SOA_COLUMN(SharedClustersNu, sharedClustersNu, uint8_t);
DECLARE_SOA_COLUMN(SharedClustersHad, sharedClustersHad, uint8_t);

DECLARE_SOA_COLUMN(IsBkgUS, isBkgUS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(ZVertex, zVertex, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, uint16_t);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);

} // namespace hadron_nuclei_tables

DECLARE_SOA_TABLE(HadronNucleiTable, "AOD", "HADNUCLEITABLE",
                  hadron_nuclei_tables::PtHad,
                  hadron_nuclei_tables::PtNu,
                  hadron_nuclei_tables::InnerParamTPCHad,
                  hadron_nuclei_tables::InnerParamTPCNu,
                  hadron_nuclei_tables::TrackIDHad,
                  hadron_nuclei_tables::TrackIDNu)
DECLARE_SOA_TABLE(HadronHyperTable, "AOD", "HADHYPERTABLE",
                  hadron_nuclei_tables::PtHyp,
                  hadron_nuclei_tables::EtaHyp,
                  hadron_nuclei_tables::PtHe3,
                  hadron_nuclei_tables::EtaHe3,
                  hadron_nuclei_tables::PhiHyp,
                  hadron_nuclei_tables::PtHad,
                  hadron_nuclei_tables::EtaHad,
                  hadron_nuclei_tables::PhiHad,
                  hadron_nuclei_tables::DcaxyHad,
                  hadron_nuclei_tables::DcazHad,
                  hadron_nuclei_tables::SignalTPCHad,
                  hadron_nuclei_tables::SignalTPCNu,
                  hadron_nuclei_tables::InnerParamTPCHad,
                  hadron_nuclei_tables::NSigmaTPCHad,
                  hadron_nuclei_tables::NSigmaTPCNu,
                  hadron_nuclei_tables::Chi2TPCHad,
                  hadron_nuclei_tables::Chi2TPCNu,
                  hadron_nuclei_tables::MassTOFHad,
                  hadron_nuclei_tables::HaddTrkHad,
                  hadron_nuclei_tables::ItsClusterSizeHad,
                  hadron_nuclei_tables::ItsClusterSizeNu,
                  hadron_nuclei_tables::SharedClustersHad,
                  hadron_nuclei_tables::TrackIDHad,
                  hadron_nuclei_tables::IsBkgUS,
                  hadron_nuclei_tables::IsBkgEM)
DECLARE_SOA_TABLE(HadronNucleiMult, "AOD", "HADNUCLEIMULT",
                  hadron_nuclei_tables::CollisionId,
                  hadron_nuclei_tables::ZVertex,
                  hadron_nuclei_tables::Multiplicity,
                  hadron_nuclei_tables::CentFT0C,
                  hadron_nuclei_tables::MultiplicityFT0C)

} // namespace o2::aod

#endif // PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_HADRONNUCLEITABLES_H_
