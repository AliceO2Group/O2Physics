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
/// \file PionNucleiTables.h
/// \brief Slim tables for piNuclei
/// \author CMY
/// \date 2025-04-10

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_PIONNUCLEITABLES_H_
#define PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_PIONNUCLEITABLES_H_

namespace o2::aod
{
namespace pion_nuclei_tables
{

DECLARE_SOA_COLUMN(PtNu, ptNu, float);
DECLARE_SOA_COLUMN(EtaNu, etaNu, float);
DECLARE_SOA_COLUMN(PhiNu, phiNu, float);
DECLARE_SOA_COLUMN(PtHyp, ptHyp, float);
DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);
DECLARE_SOA_COLUMN(EtaHyp, etaHyp, float);
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);
DECLARE_SOA_COLUMN(PhiHyp, phiHyp, float);
DECLARE_SOA_COLUMN(PtPi, ptPi, float);
DECLARE_SOA_COLUMN(EtaPi, etaPi, float);
DECLARE_SOA_COLUMN(PhiPi, phiPi, float);

DECLARE_SOA_COLUMN(DcaxyNu, dcaxyNu, float);
DECLARE_SOA_COLUMN(DcazNu, dcazNu, float);
DECLARE_SOA_COLUMN(DcaxyPi, dcaxyPi, float);
DECLARE_SOA_COLUMN(DcazPi, dcazPi, float);

DECLARE_SOA_COLUMN(SignalTPCNu, signalTPCNu, float);
DECLARE_SOA_COLUMN(InnerParamTPCNu, innerParamTPCNu, float);
DECLARE_SOA_COLUMN(SignalTPCPi, signalTPCPi, float);
DECLARE_SOA_COLUMN(InnerParamTPCPi, innerParamTPCPi, float);
DECLARE_SOA_COLUMN(NClsTPCNu, nClsTPCNu, uint8_t);
DECLARE_SOA_COLUMN(NSigmaTPCNu, nSigmaTPCNu, float);
DECLARE_SOA_COLUMN(NSigmaTPCPi, nSigmaTPCPi, float);
DECLARE_SOA_COLUMN(Chi2TPCNu, chi2TPCNu, float);
DECLARE_SOA_COLUMN(Chi2TPCPi, chi2TPCPi, float);
DECLARE_SOA_COLUMN(MassTOFNu, massTOFNu, float);
DECLARE_SOA_COLUMN(MassTOFPi, massTOFPi, float);
DECLARE_SOA_COLUMN(PidTrkNu, pidTrkNu, uint32_t);
DECLARE_SOA_COLUMN(PidTrkPi, pidTrkPi, uint32_t);
DECLARE_SOA_COLUMN(TrackIDPi, trackIDPi, int);
DECLARE_SOA_COLUMN(TrackIDNu, trackIDNu, int);

DECLARE_SOA_COLUMN(ItsClusterSizeNu, itsClusterSizeNu, uint32_t);
DECLARE_SOA_COLUMN(ItsClusterSizePi, itsClusterSizePi, uint32_t);

DECLARE_SOA_COLUMN(SharedClustersNu, sharedClustersNu, uint8_t);
DECLARE_SOA_COLUMN(SharedClustersPi, sharedClustersPi, uint8_t);

DECLARE_SOA_COLUMN(IsBkgUS, isBkgUS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(ZVertex, zVertex, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, uint16_t);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);

} // namespace pion_nuclei_tables

DECLARE_SOA_TABLE(PionNucleiTable, "AOD", "PINUCLEITABLE",
                  pion_nuclei_tables::PtPi,
                  pion_nuclei_tables::PtNu,
                  pion_nuclei_tables::InnerParamTPCPi,
                  pion_nuclei_tables::InnerParamTPCNu,
                  pion_nuclei_tables::TrackIDPi,
                  pion_nuclei_tables::TrackIDNu)
DECLARE_SOA_TABLE(PionHyperTable, "AOD", "PIHYPERTABLE",
                  pion_nuclei_tables::PtHyp,
                  pion_nuclei_tables::EtaHyp,
                  pion_nuclei_tables::PtHe3,
                  pion_nuclei_tables::EtaHe3,
                  pion_nuclei_tables::PhiHyp,
                  pion_nuclei_tables::PtPi,
                  pion_nuclei_tables::EtaPi,
                  pion_nuclei_tables::PhiPi,
                  pion_nuclei_tables::DcaxyPi,
                  pion_nuclei_tables::DcazPi,
                  pion_nuclei_tables::SignalTPCPi,
                  pion_nuclei_tables::SignalTPCNu,
                  pion_nuclei_tables::InnerParamTPCPi,
                  pion_nuclei_tables::NSigmaTPCPi,
                  pion_nuclei_tables::NSigmaTPCNu,
                  pion_nuclei_tables::Chi2TPCPi,
                  pion_nuclei_tables::Chi2TPCNu,
                  pion_nuclei_tables::MassTOFPi,
                  pion_nuclei_tables::PidTrkPi,
                  pion_nuclei_tables::ItsClusterSizePi,
                  pion_nuclei_tables::ItsClusterSizeNu,
                  pion_nuclei_tables::SharedClustersPi,
                  pion_nuclei_tables::TrackIDPi,
                  pion_nuclei_tables::IsBkgUS,
                  pion_nuclei_tables::IsBkgEM)
DECLARE_SOA_TABLE(PionNucleiMult, "AOD", "PINUCLEIMULT",
                  pion_nuclei_tables::CollisionId,
                  pion_nuclei_tables::ZVertex,
                  pion_nuclei_tables::Multiplicity,
                  pion_nuclei_tables::CentFT0C,
                  pion_nuclei_tables::MultiplicityFT0C)

} // namespace o2::aod

#endif // PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_PIONNUCLEITABLES_H_
