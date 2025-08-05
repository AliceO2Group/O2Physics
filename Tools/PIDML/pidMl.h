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

/// \file pidMl.h
/// \brief Data model for PID ML training.
///
/// \author Maja Kabus <mkabus@cern.ch>

#ifndef TOOLS_PIDML_PIDML_H_
#define TOOLS_PIDML_PIDML_H_

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace pidtracks
{
DECLARE_SOA_COLUMN(MultFV0M, multFV0M, float);                     //! Non-dynamic column with V0 multiplicity
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);                     //! Non-dynamic column with T0 multiplicity
DECLARE_SOA_COLUMN(P, p, float);                                   //! Non-dynamic column with track momentum
DECLARE_SOA_COLUMN(Px, px, float);                                 //! Non-dynamic column with track x-momentum
DECLARE_SOA_COLUMN(Py, py, float);                                 //! Non-dynamic column with track y-momentum
DECLARE_SOA_COLUMN(Pz, pz, float);                                 //! Non-dynamic column with track z-momentum
DECLARE_SOA_COLUMN(Sign, sign, float);                             //! Non-dynamic column with track sign
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, uint8_t); //!
DECLARE_SOA_COLUMN(TofExpSignalDiffEl, tofExpSignalDiffEl, float); //! Difference between signal and expected for electron
DECLARE_SOA_COLUMN(TpcExpSignalDiffEl, tpcExpSignalDiffEl, float); //! Difference between signal and expected for electron
DECLARE_SOA_COLUMN(TofExpSignalDiffMu, tofExpSignalDiffMu, float); //! Difference between signal and expected for muon
DECLARE_SOA_COLUMN(TpcExpSignalDiffMu, tpcExpSignalDiffMu, float); //! Difference between signal and expected for muon
DECLARE_SOA_COLUMN(TofExpSignalDiffPi, tofExpSignalDiffPi, float); //! Difference between signal and expected for pion
DECLARE_SOA_COLUMN(TpcExpSignalDiffPi, tpcExpSignalDiffPi, float); //! Difference between signal and expected for pion
DECLARE_SOA_COLUMN(TofExpSignalDiffKa, tofExpSignalDiffKa, float); //! Difference between signal and expected for kaon
DECLARE_SOA_COLUMN(TpcExpSignalDiffKa, tpcExpSignalDiffKa, float); //! Difference between signal and expected for kaon
DECLARE_SOA_COLUMN(TofExpSignalDiffPr, tofExpSignalDiffPr, float); //! Difference between signal and expected for proton
DECLARE_SOA_COLUMN(TpcExpSignalDiffPr, tpcExpSignalDiffPr, float); //! Difference between signal and expected for proton
} // namespace pidtracks
DECLARE_SOA_TABLE(PidTracksDataMl, "AOD", "PIDTRACKSDATAML", //! Data tracks for prediction and domain adaptation
                  aod::track::TPCSignal,
                  aod::track::TRDSignal, aod::track::TRDPattern,
                  aod::pidtofsignal::TOFSignal,
                  aod::pidtofbeta::Beta,
                  pidtracks::P,
                  aod::track::Pt,
                  pidtracks::Px,
                  pidtracks::Py,
                  pidtracks::Pz,
                  pidtracks::Sign,
                  aod::track::X,
                  aod::track::Y,
                  aod::track::Z,
                  aod::track::Alpha,
                  aod::track::TrackType,
                  aod::track::TPCNClsShared,
                  aod::track::DcaXY,
                  aod::track::DcaZ);
DECLARE_SOA_TABLE(PidTracksData, "AOD", "PIDTRACKSDATA", //! Data tracks for comparative analysis
                  aod::mult::MultFV0A, aod::mult::MultFV0C, pidtracks::MultFV0M,
                  aod::mult::MultFT0A, aod::mult::MultFT0C, pidtracks::MultFT0M,
                  aod::mult::MultZNA, aod::mult::MultZNC,
                  aod::mult::MultTracklets, aod::mult::MultTPC,
                  aod::track::TPCSignal,
                  aod::track::TRDSignal, aod::track::TRDPattern,
                  aod::track::TrackEtaEMCAL,
                  aod::track::TrackPhiEMCAL,
                  aod::pidtofsignal::TOFSignal,
                  aod::pidtofbeta::Beta,
                  pidtracks::P,
                  aod::track::Pt,
                  pidtracks::Px,
                  pidtracks::Py,
                  pidtracks::Pz,
                  pidtracks::Sign,
                  aod::track::X,
                  aod::track::Y,
                  aod::track::Z,
                  aod::track::Alpha,
                  aod::track::TrackType,
                  aod::track::TPCNClsShared,
                  aod::track::DcaXY,
                  aod::track::DcaZ,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCExpSigmaEl,
                  pidtracks::TpcExpSignalDiffEl,
                  pidtof::TOFNSigmaEl,
                  pidtof::TOFExpSigmaEl,
                  pidtracks::TofExpSignalDiffEl,
                  pidtpc::TPCNSigmaMu,
                  pidtpc::TPCExpSigmaMu,
                  pidtracks::TpcExpSignalDiffMu,
                  pidtof::TOFNSigmaMu,
                  pidtof::TOFExpSigmaMu,
                  pidtracks::TofExpSignalDiffMu,
                  pidtpc::TPCNSigmaPi,
                  pidtpc::TPCExpSigmaPi,
                  pidtracks::TpcExpSignalDiffPi,
                  pidtof::TOFNSigmaPi,
                  pidtof::TOFExpSigmaPi,
                  pidtracks::TofExpSignalDiffPi,
                  pidtpc::TPCNSigmaKa,
                  pidtpc::TPCExpSigmaKa,
                  pidtracks::TpcExpSignalDiffKa,
                  pidtof::TOFNSigmaKa,
                  pidtof::TOFExpSigmaKa,
                  pidtracks::TofExpSignalDiffKa,
                  pidtpc::TPCNSigmaPr,
                  pidtpc::TPCExpSigmaPr,
                  pidtracks::TpcExpSignalDiffPr,
                  pidtof::TOFNSigmaPr,
                  pidtof::TOFExpSigmaPr,
                  pidtracks::TofExpSignalDiffPr);
DECLARE_SOA_TABLE(PidTracksMcMl, "AOD", "PIDTRACKSMCML", //! MC tracks for training
                  aod::track::TPCSignal,
                  aod::track::TRDSignal, aod::track::TRDPattern,
                  aod::pidtofsignal::TOFSignal,
                  aod::pidtofbeta::Beta,
                  pidtracks::P,
                  aod::track::Pt,
                  pidtracks::Px,
                  pidtracks::Py,
                  pidtracks::Pz,
                  pidtracks::Sign,
                  aod::track::X,
                  aod::track::Y,
                  aod::track::Z,
                  aod::track::Alpha,
                  aod::track::TrackType,
                  aod::track::TPCNClsShared,
                  aod::track::DcaXY,
                  aod::track::DcaZ,
                  aod::mcparticle::PdgCode,
                  pidtracks::IsPhysicalPrimary);
DECLARE_SOA_TABLE(PidTracksMc, "AOD", "PIDTRACKSMC", //! MC tracks for comparative analysis
                  aod::mult::MultFV0A, aod::mult::MultFV0C, pidtracks::MultFV0M,
                  aod::mult::MultFT0A, aod::mult::MultFT0C, pidtracks::MultFT0M,
                  aod::mult::MultZNA, aod::mult::MultZNC,
                  aod::mult::MultTracklets, aod::mult::MultTPC,
                  aod::track::TPCSignal,
                  aod::track::TRDSignal, aod::track::TRDPattern,
                  aod::track::TrackEtaEMCAL,
                  aod::track::TrackPhiEMCAL,
                  aod::pidtofsignal::TOFSignal,
                  aod::pidtofbeta::Beta,
                  pidtracks::P,
                  aod::track::Pt,
                  pidtracks::Px,
                  pidtracks::Py,
                  pidtracks::Pz,
                  pidtracks::Sign,
                  aod::track::X,
                  aod::track::Y,
                  aod::track::Z,
                  aod::track::Alpha,
                  aod::track::TrackType,
                  aod::track::TPCNClsShared,
                  aod::track::DcaXY,
                  aod::track::DcaZ,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCExpSigmaEl,
                  pidtracks::TpcExpSignalDiffEl,
                  pidtof::TOFNSigmaEl,
                  pidtof::TOFExpSigmaEl,
                  pidtracks::TofExpSignalDiffEl,
                  pidtpc::TPCNSigmaMu,
                  pidtpc::TPCExpSigmaMu,
                  pidtracks::TpcExpSignalDiffMu,
                  pidtof::TOFNSigmaMu,
                  pidtof::TOFExpSigmaMu,
                  pidtracks::TofExpSignalDiffMu,
                  pidtpc::TPCNSigmaPi,
                  pidtpc::TPCExpSigmaPi,
                  pidtracks::TpcExpSignalDiffPi,
                  pidtof::TOFNSigmaPi,
                  pidtof::TOFExpSigmaPi,
                  pidtracks::TofExpSignalDiffPi,
                  pidtpc::TPCNSigmaKa,
                  pidtpc::TPCExpSigmaKa,
                  pidtracks::TpcExpSignalDiffKa,
                  pidtof::TOFNSigmaKa,
                  pidtof::TOFExpSigmaKa,
                  pidtracks::TofExpSignalDiffKa,
                  pidtpc::TPCNSigmaPr,
                  pidtpc::TPCExpSigmaPr,
                  pidtracks::TpcExpSignalDiffPr,
                  pidtof::TOFNSigmaPr,
                  pidtof::TOFExpSigmaPr,
                  pidtracks::TofExpSignalDiffPr,
                  aod::mcparticle::PdgCode,
                  pidtracks::IsPhysicalPrimary);
} // namespace o2::aod
#endif // TOOLS_PIDML_PIDML_H_
