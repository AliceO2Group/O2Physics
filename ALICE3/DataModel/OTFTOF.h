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
/// \file   OTFTOF.h
/// \author David Dobrigkeit Chinellato
/// \author Nicolo Jacazio
/// \since  11/05/2023
/// \brief  Set of tables for the ALICE3 OTFTOF information
///

#ifndef ALICE3_DATAMODEL_OTFTOF_H_
#define ALICE3_DATAMODEL_OTFTOF_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace upgrade_tof
{
DECLARE_SOA_COLUMN(NSigmaElectronInner, nSigmaElectronInner, float); //! NSigma electron InnerTOF
DECLARE_SOA_COLUMN(NSigmaMuonInner, nSigmaMuonInner, float);         //! NSigma muon InnerTOF
DECLARE_SOA_COLUMN(NSigmaPionInner, nSigmaPionInner, float);         //! NSigma pion InnerTOF
DECLARE_SOA_COLUMN(NSigmaKaonInner, nSigmaKaonInner, float);         //! NSigma kaon InnerTOF
DECLARE_SOA_COLUMN(NSigmaProtonInner, nSigmaProtonInner, float);     //! NSigma proton InnerTOF
DECLARE_SOA_COLUMN(NSigmaElectronOuter, nSigmaElectronOuter, float); //! NSigma electron OuterTOF
DECLARE_SOA_COLUMN(NSigmaMuonOuter, nSigmaMuonOuter, float);         //! NSigma muon OuterTOF
DECLARE_SOA_COLUMN(NSigmaPionOuter, nSigmaPionOuter, float);         //! NSigma pion OuterTOF
DECLARE_SOA_COLUMN(NSigmaKaonOuter, nSigmaKaonOuter, float);         //! NSigma kaon OuterTOF
DECLARE_SOA_COLUMN(NSigmaProtonOuter, nSigmaProtonOuter, float);     //! NSigma proton OuterTOF
} // namespace upgrade_tof
DECLARE_SOA_TABLE(UpgradeTofs, "AOD", "UPGRADETOF",
                  upgrade_tof::NSigmaElectronInner,
                  upgrade_tof::NSigmaMuonInner,
                  upgrade_tof::NSigmaPionInner,
                  upgrade_tof::NSigmaKaonInner,
                  upgrade_tof::NSigmaProtonInner,
                  upgrade_tof::NSigmaElectronOuter,
                  upgrade_tof::NSigmaMuonOuter,
                  upgrade_tof::NSigmaPionOuter,
                  upgrade_tof::NSigmaKaonOuter,
                  upgrade_tof::NSigmaProtonOuter);

using UpgradeTof = UpgradeTofs::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFTOF_H_
