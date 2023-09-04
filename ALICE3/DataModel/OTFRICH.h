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
/// \file   OTFRICH.h
/// \author Nicola Nicassio, University and INFN Bari
/// \since  15/05/2023
/// \brief  Set of tables for the ALICE 3 OTFRICH information
///

#ifndef ALICE3_DATAMODEL_OTFRICH_H_
#define ALICE3_DATAMODEL_OTFRICH_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace upgrade_rich
{
DECLARE_SOA_COLUMN(NSigmaElectronRich, nSigmaElectronRich, float); //! NSigma electron BarrelRich
DECLARE_SOA_COLUMN(NSigmaMuonRich, nSigmaMuonRich, float);         //! NSigma muon BarrelRich
DECLARE_SOA_COLUMN(NSigmaPionRich, nSigmaPionRich, float);         //! NSigma pion BarrelRich
DECLARE_SOA_COLUMN(NSigmaKaonRich, nSigmaKaonRich, float);         //! NSigma kaon BarrelRich
DECLARE_SOA_COLUMN(NSigmaProtonRich, nSigmaProtonRich, float);     //! NSigma proton BarrelRich
} // namespace upgrade_rich
DECLARE_SOA_TABLE(UpgradeRichs, "AOD", "UPGRADERICH",
                  upgrade_rich::NSigmaElectronRich,
                  upgrade_rich::NSigmaMuonRich,
                  upgrade_rich::NSigmaPionRich,
                  upgrade_rich::NSigmaKaonRich,
                  upgrade_rich::NSigmaProtonRich);

using UpgradeRich = UpgradeRichs::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFRICH_H_
