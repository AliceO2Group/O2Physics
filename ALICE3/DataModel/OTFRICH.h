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
DECLARE_SOA_COLUMN(NSigmaDeuteronRich, nSigmaDeuteronRich, float); //! NSigma deuteron BarrelRich
DECLARE_SOA_COLUMN(NSigmaTritonRich, nSigmaTritonRich, float);     //! NSigma triton BarrelRich
DECLARE_SOA_COLUMN(NSigmaHelium3Rich, nSigmaHelium3Rich, float);   //! NSigma helium3 BarrelRich
DECLARE_SOA_COLUMN(NSigmaAlphaRich, nSigmaAlphaRich, float);       //! NSigma alpha BarrelRich
DECLARE_SOA_DYNAMIC_COLUMN(NSigmaRich, nSigmaRich,                 //! General function to get the nSigma for the RICH
                           [](const float el,
                              const float mu,
                              const float pi,
                              const float ka,
                              const float pr,
                              const float de,
                              const float tr,
                              const float he3,
                              const float al,
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
                               case 5:
                                 return de;
                               case 6:
                                 return tr;
                               case 7:
                                 return he3;
                               case 8:
                                 return al;
                               default:
                                 LOG(fatal) << "Unrecognized PDG code for RICH";
                                 return 999.f;
                             }
                           });

DECLARE_SOA_COLUMN(HasSig, hasSig, bool);           //! Has signal in the barrel rich (is particle over threshold)
DECLARE_SOA_COLUMN(HasSigInGas, hasSigInGas, bool); //! Has signal in the gas radiator in the barrel rich (is particle over threshold)
DECLARE_SOA_COLUMN(HasSigEl, hasSigEl, bool);       //! Has nSigma electron BarrelRich (is electron over threshold)
DECLARE_SOA_COLUMN(HasSigMu, hasSigMu, bool);       //! Has nSigma muon BarrelRich (is muon over threshold)
DECLARE_SOA_COLUMN(HasSigPi, hasSigPi, bool);       //! Has nSigma pion BarrelRich (is pion over threshold)
DECLARE_SOA_COLUMN(HasSigKa, hasSigKa, bool);       //! Has nSigma kaon BarrelRich (is kaon over threshold)
DECLARE_SOA_COLUMN(HasSigPr, hasSigPr, bool);       //! Has nSigma proton BarrelRich (is proton over threshold)
DECLARE_SOA_COLUMN(HasSigDe, hasSigDe, bool);       //! Has nSigma deuteron BarrelRich (is deuteron over threshold)
DECLARE_SOA_COLUMN(HasSigTr, hasSigTr, bool);       //! Has nSigma triton BarrelRich (is triton over threshold)
DECLARE_SOA_COLUMN(HasSigHe3, hasSigHe3, bool);     //! Has nSigma helium3 BarrelRich (is helium3 over threshold)
DECLARE_SOA_COLUMN(HasSigAl, hasSigAl, bool);       //! Has nSigma alpha BarrelRich (is alpha over threshold)

} // namespace upgrade_rich
DECLARE_SOA_TABLE(UpgradeRichs, "AOD", "UPGRADERICH",
                  upgrade_rich::NSigmaElectronRich,
                  upgrade_rich::NSigmaMuonRich,
                  upgrade_rich::NSigmaPionRich,
                  upgrade_rich::NSigmaKaonRich,
                  upgrade_rich::NSigmaProtonRich,
                  upgrade_rich::NSigmaDeuteronRich,
                  upgrade_rich::NSigmaTritonRich,
                  upgrade_rich::NSigmaHelium3Rich,
                  upgrade_rich::NSigmaAlphaRich,
                  upgrade_rich::NSigmaRich<upgrade_rich::NSigmaElectronRich,
                                           upgrade_rich::NSigmaMuonRich,
                                           upgrade_rich::NSigmaPionRich,
                                           upgrade_rich::NSigmaKaonRich,
                                           upgrade_rich::NSigmaProtonRich,
                                           upgrade_rich::NSigmaDeuteronRich,
                                           upgrade_rich::NSigmaTritonRich,
                                           upgrade_rich::NSigmaHelium3Rich,
                                           upgrade_rich::NSigmaAlphaRich>);

using UpgradeRich = UpgradeRichs::iterator;

DECLARE_SOA_TABLE(UpgradeRichSignals, "AOD", "UPGRADERICHSIG",
                  upgrade_rich::HasSig,
                  upgrade_rich::HasSigEl,
                  upgrade_rich::HasSigMu,
                  upgrade_rich::HasSigPi,
                  upgrade_rich::HasSigKa,
                  upgrade_rich::HasSigPr,
                  upgrade_rich::HasSigDe,
                  upgrade_rich::HasSigTr,
                  upgrade_rich::HasSigHe3,
                  upgrade_rich::HasSigAl,
                  upgrade_rich::HasSigInGas);

using UpgradeRichSignal = UpgradeRichSignals::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFRICH_H_
