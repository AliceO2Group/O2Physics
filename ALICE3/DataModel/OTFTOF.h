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
DECLARE_SOA_COLUMN(NSigmaElectronInnerTOF, nSigmaElectronInnerTOF, float);     //! NSigma electron InnerTOF
DECLARE_SOA_COLUMN(NSigmaMuonInnerTOF, nSigmaMuonInnerTOF, float);             //! NSigma muon InnerTOF
DECLARE_SOA_COLUMN(NSigmaPionInnerTOF, nSigmaPionInnerTOF, float);             //! NSigma pion InnerTOF
DECLARE_SOA_COLUMN(NSigmaKaonInnerTOF, nSigmaKaonInnerTOF, float);             //! NSigma kaon InnerTOF
DECLARE_SOA_COLUMN(NSigmaProtonInnerTOF, nSigmaProtonInnerTOF, float);         //! NSigma proton InnerTOF
DECLARE_SOA_COLUMN(InnerTOFTrackLength, innerTOFTrackLength, float);           //! track length for calculation of InnerTOF
DECLARE_SOA_COLUMN(InnerTOFTrackLengthReco, innerTOFTrackLengthReco, float);   //! track length for calculation of InnerTOF
DECLARE_SOA_COLUMN(DeltaTrackLengthInnerTOF, deltaTrackLengthInnerTOF, float); //! track length for calculation of InnerTOF
DECLARE_SOA_COLUMN(NSigmaElectronOuterTOF, nSigmaElectronOuterTOF, float);     //! NSigma electron OuterTOF
DECLARE_SOA_COLUMN(NSigmaMuonOuterTOF, nSigmaMuonOuterTOF, float);             //! NSigma muon OuterTOF
DECLARE_SOA_COLUMN(NSigmaPionOuterTOF, nSigmaPionOuterTOF, float);             //! NSigma pion OuterTOF
DECLARE_SOA_COLUMN(NSigmaKaonOuterTOF, nSigmaKaonOuterTOF, float);             //! NSigma kaon OuterTOF
DECLARE_SOA_COLUMN(NSigmaProtonOuterTOF, nSigmaProtonOuterTOF, float);         //! NSigma proton OuterTOF
DECLARE_SOA_COLUMN(OuterTOFTrackLength, outerTOFTrackLength, float);           //! track length for calculation of OuterTOF
DECLARE_SOA_COLUMN(OuterTOFTrackLengthReco, outerTOFTrackLengthReco, float);   //! track length for calculation of OuterTOF
DECLARE_SOA_COLUMN(DeltaTrackLengthOuterTOF, deltaTrackLengthOuterTOF, float); //! track length for calculation of InnerTOF
} // namespace upgrade_tof
DECLARE_SOA_TABLE(UpgradeTofs, "AOD", "UPGRADETOF",
                  upgrade_tof::NSigmaElectronInnerTOF,
                  upgrade_tof::NSigmaMuonInnerTOF,
                  upgrade_tof::NSigmaPionInnerTOF,
                  upgrade_tof::NSigmaKaonInnerTOF,
                  upgrade_tof::NSigmaProtonInnerTOF,
                  upgrade_tof::InnerTOFTrackLength,
                  upgrade_tof::InnerTOFTrackLengthReco,
                  upgrade_tof::DeltaTrackLengthInnerTOF,
                  upgrade_tof::NSigmaElectronOuterTOF,
                  upgrade_tof::NSigmaMuonOuterTOF,
                  upgrade_tof::NSigmaPionOuterTOF,
                  upgrade_tof::NSigmaKaonOuterTOF,
                  upgrade_tof::NSigmaProtonOuterTOF,
                  upgrade_tof::OuterTOFTrackLength,
                  upgrade_tof::OuterTOFTrackLengthReco,
                  upgrade_tof::DeltaTrackLengthOuterTOF);

using UpgradeTof = UpgradeTofs::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFTOF_H_
