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
/// \file   LFParticleIdentification.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022/11/20
/// \brief  Set of tables for the particle identification, LF oriented.
///

#ifndef PWGLF_DATAMODEL_LFPARTICLEIDENTIFICATION_H_
#define PWGLF_DATAMODEL_LFPARTICLEIDENTIFICATION_H_

#include "Common/DataModel/PIDResponseTPC.h"

namespace o2::aod
{
DECLARE_SOA_TABLE(pidTPCLfFullEl, "AOD", "pidTPCLfFullEl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for electron
                  o2::soa::Marker<1>,
                  pidtpc::TPCExpSignalEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>, pidtpc::TPCExpSignalDiffEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>, pidtpc::TPCExpSigmaEl, pidtpc::TPCNSigmaEl);
DECLARE_SOA_TABLE(pidTPCLfFullMu, "AOD", "pidTPCLfFullMu", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for muon
                  o2::soa::Marker<2>,
                  pidtpc::TPCExpSignalMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>, pidtpc::TPCExpSignalDiffMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>, pidtpc::TPCExpSigmaMu, pidtpc::TPCNSigmaMu);
DECLARE_SOA_TABLE(pidTPCLfFullPi, "AOD", "pidTPCLfFullPi", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for pion
                  o2::soa::Marker<3>,
                  pidtpc::TPCExpSignalPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>, pidtpc::TPCExpSignalDiffPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>, pidtpc::TPCExpSigmaPi, pidtpc::TPCNSigmaPi);
DECLARE_SOA_TABLE(pidTPCLfFullKa, "AOD", "pidTPCLfFullKa", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for kaon
                  o2::soa::Marker<4>,
                  pidtpc::TPCExpSignalKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>, pidtpc::TPCExpSignalDiffKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>, pidtpc::TPCExpSigmaKa, pidtpc::TPCNSigmaKa);
DECLARE_SOA_TABLE(pidTPCLfFullPr, "AOD", "pidTPCLfFullPr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for proton
                  o2::soa::Marker<5>,
                  pidtpc::TPCExpSignalPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>, pidtpc::TPCExpSignalDiffPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>, pidtpc::TPCExpSigmaPr, pidtpc::TPCNSigmaPr);
DECLARE_SOA_TABLE(pidTPCLfFullDe, "AOD", "pidTPCLfFullDe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for deuteron
                  o2::soa::Marker<6>,
                  pidtpc::TPCExpSignalDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>, pidtpc::TPCExpSignalDiffDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>, pidtpc::TPCExpSigmaDe, pidtpc::TPCNSigmaDe);
DECLARE_SOA_TABLE(pidTPCLfFullTr, "AOD", "pidTPCLfFullTr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for triton
                  o2::soa::Marker<7>,
                  pidtpc::TPCExpSignalTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>, pidtpc::TPCExpSignalDiffTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>, pidtpc::TPCExpSigmaTr, pidtpc::TPCNSigmaTr);
DECLARE_SOA_TABLE(pidTPCLfFullHe, "AOD", "pidTPCLfFullHe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for helium3
                  o2::soa::Marker<8>,
                  pidtpc::TPCExpSignalHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>, pidtpc::TPCExpSignalDiffHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>, pidtpc::TPCExpSigmaHe, pidtpc::TPCNSigmaHe);
DECLARE_SOA_TABLE(pidTPCLfFullAl, "AOD", "pidTPCLfFullAl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for alpha
                  o2::soa::Marker<9>,
                  pidtpc::TPCExpSignalAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>, pidtpc::TPCExpSignalDiffAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>, pidtpc::TPCExpSigmaAl, pidtpc::TPCNSigmaAl);

// Tiny size tables
DECLARE_SOA_TABLE(pidTPCLfEl, "AOD", "pidTPCLfEl", //! Table of the TPC response with binned Nsigma for electron
                  o2::soa::Marker<10>,
                  pidtpc_tiny::TPCNSigmaStoreEl, pidtpc_tiny::TPCNSigmaEl<pidtpc_tiny::TPCNSigmaStoreEl>);
DECLARE_SOA_TABLE(pidTPCLfMu, "AOD", "pidTPCLfMu", //! Table of the TPC response with binned Nsigma for muon
                  o2::soa::Marker<11>,
                  pidtpc_tiny::TPCNSigmaStoreMu, pidtpc_tiny::TPCNSigmaMu<pidtpc_tiny::TPCNSigmaStoreMu>);
DECLARE_SOA_TABLE(pidTPCLfPi, "AOD", "pidTPCLfPi", //! Table of the TPC response with binned Nsigma for pion
                  o2::soa::Marker<12>,
                  pidtpc_tiny::TPCNSigmaStorePi, pidtpc_tiny::TPCNSigmaPi<pidtpc_tiny::TPCNSigmaStorePi>);
DECLARE_SOA_TABLE(pidTPCLfKa, "AOD", "pidTPCLfKa", //! Table of the TPC response with binned Nsigma for kaon
                  o2::soa::Marker<13>,
                  pidtpc_tiny::TPCNSigmaStoreKa, pidtpc_tiny::TPCNSigmaKa<pidtpc_tiny::TPCNSigmaStoreKa>);
DECLARE_SOA_TABLE(pidTPCLfPr, "AOD", "pidTPCLfPr", //! Table of the TPC response with binned Nsigma for proton
                  o2::soa::Marker<14>,
                  pidtpc_tiny::TPCNSigmaStorePr, pidtpc_tiny::TPCNSigmaPr<pidtpc_tiny::TPCNSigmaStorePr>);
DECLARE_SOA_TABLE(pidTPCLfDe, "AOD", "pidTPCLfDe", //! Table of the TPC response with binned Nsigma for deuteron
                  o2::soa::Marker<15>,
                  pidtpc_tiny::TPCNSigmaStoreDe, pidtpc_tiny::TPCNSigmaDe<pidtpc_tiny::TPCNSigmaStoreDe>);
DECLARE_SOA_TABLE(pidTPCLfTr, "AOD", "pidTPCLfTr", //! Table of the TPC response with binned Nsigma for triton
                  o2::soa::Marker<16>,
                  pidtpc_tiny::TPCNSigmaStoreTr, pidtpc_tiny::TPCNSigmaTr<pidtpc_tiny::TPCNSigmaStoreTr>);
DECLARE_SOA_TABLE(pidTPCLfHe, "AOD", "pidTPCLfHe", //! Table of the TPC response with binned Nsigma for helium3
                  o2::soa::Marker<17>,
                  pidtpc_tiny::TPCNSigmaStoreHe, pidtpc_tiny::TPCNSigmaHe<pidtpc_tiny::TPCNSigmaStoreHe>);
DECLARE_SOA_TABLE(pidTPCLfAl, "AOD", "pidTPCLfAl", //! Table of the TPC response with binned Nsigma for alpha
                  o2::soa::Marker<18>,
                  pidtpc_tiny::TPCNSigmaStoreAl, pidtpc_tiny::TPCNSigmaAl<pidtpc_tiny::TPCNSigmaStoreAl>);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFPARTICLEIDENTIFICATION_H_
