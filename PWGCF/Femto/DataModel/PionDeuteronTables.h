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
/// \file PionDeuteronTables.h
/// \brief Slim tables for piDeuteron
/// \author CMY
/// \date 2025-04-10

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGCF_FEMTO_DATAMODEL_PIONDEUTERONTABLES_H_
#define PWGCF_FEMTO_DATAMODEL_PIONDEUTERONTABLES_H_

namespace o2::aod
{
namespace pion_deuteron_tables
{

DECLARE_SOA_COLUMN(PtDe, ptDe, float);
DECLARE_SOA_COLUMN(EtaDe, etaDe, float);
DECLARE_SOA_COLUMN(PhiDe, phiDe, float);
DECLARE_SOA_COLUMN(PtPi, ptPi, float);
DECLARE_SOA_COLUMN(EtaPi, etaPi, float);
DECLARE_SOA_COLUMN(PhiPi, phiPi, float);

DECLARE_SOA_COLUMN(DcaxyDe, dcaxyDe, float);
DECLARE_SOA_COLUMN(DcazDe, dcazDe, float);
DECLARE_SOA_COLUMN(DcaxyPi, dcaxyPi, float);
DECLARE_SOA_COLUMN(DcazPi, dcazPi, float);

DECLARE_SOA_COLUMN(SignalTPCDe, signalTPCDe, float);
DECLARE_SOA_COLUMN(InnerParamTPCDe, innerParamTPCDe, float);
DECLARE_SOA_COLUMN(SignalTPCPi, signalTPCPi, float);
DECLARE_SOA_COLUMN(InnerParamTPCPi, innerParamTPCPi, float);
DECLARE_SOA_COLUMN(NClsTPCDe, nClsTPCDe, uint8_t);
DECLARE_SOA_COLUMN(NSigmaTPCDe, nSigmaTPCDe, float);
DECLARE_SOA_COLUMN(NSigmaTPCPi, nSigmaTPCPi, float);
DECLARE_SOA_COLUMN(Chi2TPCDe, chi2TPCDe, float);
DECLARE_SOA_COLUMN(Chi2TPCPi, chi2TPCPi, float);
DECLARE_SOA_COLUMN(MassTOFDe, massTOFDe, float);
DECLARE_SOA_COLUMN(MassTOFPi, massTOFPi, float);
DECLARE_SOA_COLUMN(PidTrkDe, pidTrkDe, uint32_t);
DECLARE_SOA_COLUMN(PidTrkPi, pidTrkPi, uint32_t);

DECLARE_SOA_COLUMN(ItsClusterSizeDe, itsClusterSizeDe, uint32_t);
DECLARE_SOA_COLUMN(ItsClusterSizePi, itsClusterSizePi, uint32_t);

DECLARE_SOA_COLUMN(SharedClustersDe, sharedClustersDe, uint8_t);
DECLARE_SOA_COLUMN(SharedClustersPi, sharedClustersPi, uint8_t);

DECLARE_SOA_COLUMN(IsBkgUS, isBkgUS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(ZVertex, zVertex, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, uint16_t);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);

} // namespace pion_deuteron_tables

DECLARE_SOA_TABLE(PionDeuteronTable, "AOD", "PIDEUTERONTABLE",
                  pion_deuteron_tables::PtDe,
                  pion_deuteron_tables::EtaDe,
                  pion_deuteron_tables::PhiDe,
                  pion_deuteron_tables::PtPi,
                  pion_deuteron_tables::EtaPi,
                  pion_deuteron_tables::PhiPi,
                  pion_deuteron_tables::DcaxyDe,
                  pion_deuteron_tables::DcazDe,
                  pion_deuteron_tables::DcaxyPi,
                  pion_deuteron_tables::DcazPi,
                  pion_deuteron_tables::SignalTPCDe,
                  pion_deuteron_tables::InnerParamTPCDe,
                  pion_deuteron_tables::SignalTPCPi,
                  pion_deuteron_tables::InnerParamTPCPi,
                  pion_deuteron_tables::NClsTPCDe,
                  pion_deuteron_tables::NSigmaTPCDe,
                  pion_deuteron_tables::NSigmaTPCPi,
                  pion_deuteron_tables::Chi2TPCDe,
                  pion_deuteron_tables::Chi2TPCPi,
                  pion_deuteron_tables::MassTOFDe,
                  pion_deuteron_tables::MassTOFPi,
                  pion_deuteron_tables::PidTrkDe,
                  pion_deuteron_tables::PidTrkPi,
                  pion_deuteron_tables::ItsClusterSizeDe,
                  pion_deuteron_tables::ItsClusterSizePi,
                  pion_deuteron_tables::SharedClustersDe,
                  pion_deuteron_tables::SharedClustersPi,
                  pion_deuteron_tables::IsBkgUS,
                  pion_deuteron_tables::IsBkgEM)
DECLARE_SOA_TABLE(PionDeuteronMult, "AOD", "PIDEUTERONMULT",
                  pion_deuteron_tables::CollisionId,
                  pion_deuteron_tables::ZVertex,
                  pion_deuteron_tables::Multiplicity,
                  pion_deuteron_tables::CentFT0C,
                  pion_deuteron_tables::MultiplicityFT0C)

} // namespace o2::aod

#endif // PWGCF_FEMTO_DATAMODEL_PIONDEUTERONTABLES_H_
