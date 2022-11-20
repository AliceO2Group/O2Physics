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

// O2 includes
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "ReconstructionDataFormats/PID.h"
#include "Framework/Logger.h"

namespace o2::aod
{
// DECLARE_SOA_TABLE(lfpidTPCFullEl, "AOD", "lfpidTPCFullEl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for electron
//                   pidtpc::TPCExpSignalEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>, pidtpc::TPCExpSignalDiffEl<pidtpc::TPCNSigmaEl, pidtpc::TPCExpSigmaEl>, pidtpc::TPCExpSigmaEl, pidtpc::TPCNSigmaEl);
// DECLARE_SOA_TABLE(lfpidTPCFullMu, "AOD", "lfpidTPCFullMu", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for muon
//                   pidtpc::TPCExpSignalMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>, pidtpc::TPCExpSignalDiffMu<pidtpc::TPCNSigmaMu, pidtpc::TPCExpSigmaMu>, pidtpc::TPCExpSigmaMu, pidtpc::TPCNSigmaMu);
// DECLARE_SOA_TABLE(lfpidTPCFullPi, "AOD", "lfpidTPCFullPi", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for pion
//                   pidtpc::TPCExpSignalPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>, pidtpc::TPCExpSignalDiffPi<pidtpc::TPCNSigmaPi, pidtpc::TPCExpSigmaPi>, pidtpc::TPCExpSigmaPi, pidtpc::TPCNSigmaPi);
// DECLARE_SOA_TABLE(lfpidTPCFullKa, "AOD", "lfpidTPCFullKa", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for kaon
//                   pidtpc::TPCExpSignalKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>, pidtpc::TPCExpSignalDiffKa<pidtpc::TPCNSigmaKa, pidtpc::TPCExpSigmaKa>, pidtpc::TPCExpSigmaKa, pidtpc::TPCNSigmaKa);
// DECLARE_SOA_TABLE(lfpidTPCFullPr, "AOD", "lfpidTPCFullPr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for proton
//                   pidtpc::TPCExpSignalPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>, pidtpc::TPCExpSignalDiffPr<pidtpc::TPCNSigmaPr, pidtpc::TPCExpSigmaPr>, pidtpc::TPCExpSigmaPr, pidtpc::TPCNSigmaPr);
// DECLARE_SOA_TABLE(lfpidTPCFullDe, "AOD", "lfpidTPCFullDe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for deuteron
//                   pidtpc::TPCExpSignalDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>, pidtpc::TPCExpSignalDiffDe<pidtpc::TPCNSigmaDe, pidtpc::TPCExpSigmaDe>, pidtpc::TPCExpSigmaDe, pidtpc::TPCNSigmaDe);
// DECLARE_SOA_TABLE(lfpidTPCFullTr, "AOD", "lfpidTPCFullTr", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for triton
//                   pidtpc::TPCExpSignalTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>, pidtpc::TPCExpSignalDiffTr<pidtpc::TPCNSigmaTr, pidtpc::TPCExpSigmaTr>, pidtpc::TPCExpSigmaTr, pidtpc::TPCNSigmaTr);
// DECLARE_SOA_TABLE(lfpidTPCFullHe, "AOD", "lfpidTPCFullHe", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for helium3
//                   pidtpc::TPCExpSignalHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>, pidtpc::TPCExpSignalDiffHe<pidtpc::TPCNSigmaHe, pidtpc::TPCExpSigmaHe>, pidtpc::TPCExpSigmaHe, pidtpc::TPCNSigmaHe);
// DECLARE_SOA_TABLE(lfpidTPCFullAl, "AOD", "lfpidTPCFullAl", //! Table of the TPC (full) response with expected signal, expected resolution and Nsigma for alpha
//                   pidtpc::TPCExpSignalAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>, pidtpc::TPCExpSignalDiffAl<pidtpc::TPCNSigmaAl, pidtpc::TPCExpSigmaAl>, pidtpc::TPCExpSigmaAl, pidtpc::TPCNSigmaAl);

// // Tiny size tables
// DECLARE_SOA_TABLE(lfpidTPCEl, "AOD", "lfpidTPCEl", //! Table of the TPC response with binned Nsigma for electron
//                   pidtpc_tiny::TPCNSigmaStoreEl, pidtpc_tiny::TPCNSigmaEl<pidtpc_tiny::TPCNSigmaStoreEl>);
// DECLARE_SOA_TABLE(lfpidTPCMu, "AOD", "lfpidTPCMu", //! Table of the TPC response with binned Nsigma for muon
//                   pidtpc_tiny::TPCNSigmaStoreMu, pidtpc_tiny::TPCNSigmaMu<pidtpc_tiny::TPCNSigmaStoreMu>);
// DECLARE_SOA_TABLE(lfpidTPCPi, "AOD", "lfpidTPCPi", //! Table of the TPC response with binned Nsigma for pion
//                   pidtpc_tiny::TPCNSigmaStorePi, pidtpc_tiny::TPCNSigmaPi<pidtpc_tiny::TPCNSigmaStorePi>);
// DECLARE_SOA_TABLE(lfpidTPCKa, "AOD", "lfpidTPCKa", //! Table of the TPC response with binned Nsigma for kaon
//                   pidtpc_tiny::TPCNSigmaStoreKa, pidtpc_tiny::TPCNSigmaKa<pidtpc_tiny::TPCNSigmaStoreKa>);
// DECLARE_SOA_TABLE(lfpidTPCPr, "AOD", "lfpidTPCPr", //! Table of the TPC response with binned Nsigma for proton
//                   pidtpc_tiny::TPCNSigmaStorePr, pidtpc_tiny::TPCNSigmaPr<pidtpc_tiny::TPCNSigmaStorePr>);
// DECLARE_SOA_TABLE(lfpidTPCDe, "AOD", "lfpidTPCDe", //! Table of the TPC response with binned Nsigma for deuteron
//                   pidtpc_tiny::TPCNSigmaStoreDe, pidtpc_tiny::TPCNSigmaDe<pidtpc_tiny::TPCNSigmaStoreDe>);
// DECLARE_SOA_TABLE(lfpidTPCTr, "AOD", "lfpidTPCTr", //! Table of the TPC response with binned Nsigma for triton
//                   pidtpc_tiny::TPCNSigmaStoreTr, pidtpc_tiny::TPCNSigmaTr<pidtpc_tiny::TPCNSigmaStoreTr>);
// DECLARE_SOA_TABLE(lfpidTPCHe, "AOD", "lfpidTPCHe", //! Table of the TPC response with binned Nsigma for helium3
//                   pidtpc_tiny::TPCNSigmaStoreHe, pidtpc_tiny::TPCNSigmaHe<pidtpc_tiny::TPCNSigmaStoreHe>);
// DECLARE_SOA_TABLE(lfpidTPCAl, "AOD", "lfpidTPCAl", //! Table of the TPC response with binned Nsigma for alpha
//                   pidtpc_tiny::TPCNSigmaStoreAl, pidtpc_tiny::TPCNSigmaAl<pidtpc_tiny::TPCNSigmaStoreAl>);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFPARTICLEIDENTIFICATION_H_
