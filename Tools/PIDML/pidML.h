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
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/PID/PIDResponse.h"

using namespace o2;

namespace o2::aod
{
namespace pidtracks
{
DECLARE_SOA_COLUMN(MultV0M, multV0M, float);                       //! Non-dynamic column with V0 multiplicity
DECLARE_SOA_COLUMN(MultT0M, multT0M, float);                       //! Non-dynamic column with T0 multiplicity
DECLARE_SOA_COLUMN(P, p, float);                                   //! Non-dynamic column with track momentum
DECLARE_SOA_COLUMN(Px, px, float);                                 //! Non-dynamic column with track x-momentum
DECLARE_SOA_COLUMN(Py, py, float);                                 //! Non-dynamic column with track y-momentum
DECLARE_SOA_COLUMN(Pz, pz, float);                                 //! Non-dynamic column with track z-momentum
DECLARE_SOA_COLUMN(Sign, sign, float);                             //! Non-dynamic column with track sign
DECLARE_SOA_COLUMN(TOFSignal, tofSignal, float);                   //! Private version of the TOF signal
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, uint8_t); //!
DECLARE_SOA_COLUMN(TOFExpSignalDiffEl, tofExpSignalDiffEl, float); //! Difference between signal and expected for electron
DECLARE_SOA_COLUMN(TPCExpSignalDiffEl, tpcExpSignalDiffEl, float); //! Difference between signal and expected for electron
DECLARE_SOA_COLUMN(TOFExpSignalDiffMu, tofExpSignalDiffMu, float); //! Difference between signal and expected for muon
DECLARE_SOA_COLUMN(TPCExpSignalDiffMu, tpcExpSignalDiffMu, float); //! Difference between signal and expected for muon
DECLARE_SOA_COLUMN(TOFExpSignalDiffPi, tofExpSignalDiffPi, float); //! Difference between signal and expected for pion
DECLARE_SOA_COLUMN(TPCExpSignalDiffPi, tpcExpSignalDiffPi, float); //! Difference between signal and expected for pion
DECLARE_SOA_COLUMN(TOFExpSignalDiffKa, tofExpSignalDiffKa, float); //! Difference between signal and expected for kaon
DECLARE_SOA_COLUMN(TPCExpSignalDiffKa, tpcExpSignalDiffKa, float); //! Difference between signal and expected for kaon
DECLARE_SOA_COLUMN(TOFExpSignalDiffPr, tofExpSignalDiffPr, float); //! Difference between signal and expected for proton
DECLARE_SOA_COLUMN(TPCExpSignalDiffPr, tpcExpSignalDiffPr, float); //! Difference between signal and expected for proton
} // namespace pidtracks
DECLARE_SOA_TABLE(PidTracksReal, "AOD", "PIDTRACKSREAL", //! Real tracks for prediction and domain adaptation
                  aod::cent::CentEstV0M,
                  aod::mult::MultV0A, aod::mult::MultV0C, pidtracks::MultV0M,
                  aod::mult::MultT0A, aod::mult::MultT0C, pidtracks::MultT0M,
                  aod::mult::MultZNA, aod::mult::MultZNC,
                  aod::mult::MultTracklets, aod::mult::MultTPC,
                  aod::track::TPCSignal,
                  aod::track::TRDSignal,
                  aod::track::TrackEtaEMCAL,
                  aod::track::TrackPhiEMCAL,
                  aod::pidtracks::TOFSignal,
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
                  pidtracks::TPCExpSignalDiffEl,
                  pidtof::TOFNSigmaEl,
                  pidtof::TOFExpSigmaEl,
                  pidtracks::TOFExpSignalDiffEl,
                  pidtpc::TPCNSigmaMu,
                  pidtpc::TPCExpSigmaMu,
                  pidtracks::TPCExpSignalDiffMu,
                  pidtof::TOFNSigmaMu,
                  pidtof::TOFExpSigmaMu,
                  pidtracks::TOFExpSignalDiffMu,
                  pidtpc::TPCNSigmaPi,
                  pidtpc::TPCExpSigmaPi,
                  pidtracks::TPCExpSignalDiffPi,
                  pidtof::TOFNSigmaPi,
                  pidtof::TOFExpSigmaPi,
                  pidtracks::TOFExpSignalDiffPi,
                  pidtpc::TPCNSigmaKa,
                  pidtpc::TPCExpSigmaKa,
                  pidtracks::TPCExpSignalDiffKa,
                  pidtof::TOFNSigmaKa,
                  pidtof::TOFExpSigmaKa,
                  pidtracks::TOFExpSignalDiffKa,
                  pidtpc::TPCNSigmaPr,
                  pidtpc::TPCExpSigmaPr,
                  pidtracks::TPCExpSignalDiffPr,
                  pidtof::TOFNSigmaPr,
                  pidtof::TOFExpSigmaPr,
                  pidtracks::TOFExpSignalDiffPr);
DECLARE_SOA_TABLE(PidTracksMc, "AOD", "PIDTRACKSMC", //! MC tracks for training
                  aod::cent::CentEstV0M,
                  aod::mult::MultV0A, aod::mult::MultV0C, pidtracks::MultV0M,
                  aod::mult::MultT0A, aod::mult::MultT0C, pidtracks::MultT0M,
                  aod::mult::MultZNA, aod::mult::MultZNC,
                  aod::mult::MultTracklets, aod::mult::MultTPC,
                  aod::track::TPCSignal,
                  aod::track::TRDSignal,
                  aod::track::TrackEtaEMCAL,
                  aod::track::TrackPhiEMCAL,
                  aod::pidtracks::TOFSignal,
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
                  pidtracks::TPCExpSignalDiffEl,
                  pidtof::TOFNSigmaEl,
                  pidtof::TOFExpSigmaEl,
                  pidtracks::TOFExpSignalDiffEl,
                  pidtpc::TPCNSigmaMu,
                  pidtpc::TPCExpSigmaMu,
                  pidtracks::TPCExpSignalDiffMu,
                  pidtof::TOFNSigmaMu,
                  pidtof::TOFExpSigmaMu,
                  pidtracks::TOFExpSignalDiffMu,
                  pidtpc::TPCNSigmaPi,
                  pidtpc::TPCExpSigmaPi,
                  pidtracks::TPCExpSignalDiffPi,
                  pidtof::TOFNSigmaPi,
                  pidtof::TOFExpSigmaPi,
                  pidtracks::TOFExpSignalDiffPi,
                  pidtpc::TPCNSigmaKa,
                  pidtpc::TPCExpSigmaKa,
                  pidtracks::TPCExpSignalDiffKa,
                  pidtof::TOFNSigmaKa,
                  pidtof::TOFExpSigmaKa,
                  pidtracks::TOFExpSignalDiffKa,
                  pidtpc::TPCNSigmaPr,
                  pidtpc::TPCExpSigmaPr,
                  pidtracks::TPCExpSignalDiffPr,
                  pidtof::TOFNSigmaPr,
                  pidtof::TOFExpSigmaPr,
                  pidtracks::TOFExpSignalDiffPr,
                  aod::mcparticle::PdgCode,
                  pidtracks::IsPhysicalPrimary);
} // namespace o2::aod
