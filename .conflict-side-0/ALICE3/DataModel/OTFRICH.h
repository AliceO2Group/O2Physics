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
DECLARE_SOA_DYNAMIC_COLUMN(NSigmaRich, nSigmaRich,                 //! General function to get the nSigma for the RICH
                           [](const float el,
                              const float mu,
                              const float pi,
                              const float ka,
                              const float pr,
                              const int id) -> float {
                             switch (std::abs(id)) {
                               case 0:
                                 return el;
                               case 1:
                                 return mu;
                               case 2:
                                 return pi;
                               case 3:
                                 return ka;
                               case 4:
                                 return pr;
                               default:
                                 LOG(fatal) << "Unrecognized PDG code for RICH";
                                 return 999.f;
                             }
                           });

DECLARE_SOA_COLUMN(HasSig, hasSig, bool);     //! Has signal in the barrel rich (is particle over threshold)
DECLARE_SOA_COLUMN(HasSigEl, hasSigEl, bool); //! Has nSigma electron BarrelRich (is electron over threshold)
DECLARE_SOA_COLUMN(HasSigMu, hasSigMu, bool); //! Has nSigma muon BarrelRich (is muon over threshold)
DECLARE_SOA_COLUMN(HasSigPi, hasSigPi, bool); //! Has nSigma pion BarrelRich (is pion over threshold)
DECLARE_SOA_COLUMN(HasSigKa, hasSigKa, bool); //! Has nSigma kaon BarrelRich (is kaon over threshold)
DECLARE_SOA_COLUMN(HasSigPr, hasSigPr, bool); //! Has nSigma proton BarrelRich (is proton over threshold)

} // namespace upgrade_rich
DECLARE_SOA_TABLE(UpgradeRichs, "AOD", "UPGRADERICH",
                  upgrade_rich::NSigmaElectronRich,
                  upgrade_rich::NSigmaMuonRich,
                  upgrade_rich::NSigmaPionRich,
                  upgrade_rich::NSigmaKaonRich,
                  upgrade_rich::NSigmaProtonRich,
                  upgrade_rich::NSigmaRich<upgrade_rich::NSigmaElectronRich,
                                           upgrade_rich::NSigmaMuonRich,
                                           upgrade_rich::NSigmaPionRich,
                                           upgrade_rich::NSigmaKaonRich,
                                           upgrade_rich::NSigmaProtonRich>);

using UpgradeRich = UpgradeRichs::iterator;

DECLARE_SOA_TABLE(UpgradeRichSignals, "AOD", "UPGRADERICHSIG",
                  upgrade_rich::HasSig,
                  upgrade_rich::HasSigEl,
                  upgrade_rich::HasSigMu,
                  upgrade_rich::HasSigPi,
                  upgrade_rich::HasSigKa,
                  upgrade_rich::HasSigPr);

using UpgradeRichSignal = UpgradeRichSignals::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFRICH_H_
