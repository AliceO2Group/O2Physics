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
/// \file   tofSkimsTableCreator.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  19/11/2022
/// \brief  Definition of the skimmed data format for the TOF skims
///

#ifndef DPG_TASKS_TOF_TOFSKIMSTABLECREATOR_H_
#define DPG_TASKS_TOF_TOFSKIMSTABLECREATOR_H_

#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace tofskims
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);              //! Index to the collision
DECLARE_SOA_COLUMN(P, p, float);                             //! Momentum of the track
DECLARE_SOA_COLUMN(PtSigned, ptSigned, float);               //! Pt (signed) of the track
DECLARE_SOA_COLUMN(DeltaP, deltaP, float);                   //! Momentum difference with respect to the reference track
DECLARE_SOA_COLUMN(DeltaPt, deltaPt, float);                 //! Transverse momentum difference with respect to the reference track
DECLARE_SOA_COLUMN(Eta, eta, float);                         //! Eta of the track
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);               //! Eta difference with respect to the reference track
DECLARE_SOA_COLUMN(Phi, phi, float);                         //! Phi of the track
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);               //! Phi difference with respect to the reference track
DECLARE_SOA_COLUMN(DeltaTEl, deltaTEl, float);               //! Difference between the tofSignal and the electron expected signal
DECLARE_SOA_COLUMN(DeltaTMu, deltaTMu, float);               //! Difference between the tofSignal and the muon expected signal
DECLARE_SOA_COLUMN(DeltaTPi, deltaTPi, float);               //! Difference between the tofSignal and the pion expected signal
DECLARE_SOA_COLUMN(DeltaTKa, deltaTKa, float);               //! Difference between the tofSignal and the kaon expected signal
DECLARE_SOA_COLUMN(DeltaTPr, deltaTPr, float);               //! Difference between the tofSignal and the proton expected signal
DECLARE_SOA_COLUMN(DoubleDelta, doubleDelta, float);         //! Double difference between DeltaT
DECLARE_SOA_COLUMN(RefSign, refSign, short);                 //! Sign of the reference track
DECLARE_SOA_COLUMN(PIDForTracking, pidForTracking, uint8_t); //! Index for mass hypothesis used in tracking see PID.h for definition
DECLARE_SOA_COLUMN(EvTimeT0AC, evTimeT0AC, float);           //! Event time of the track computed with the T0AC
DECLARE_SOA_COLUMN(EvTimeT0ACErr, evTimeT0ACErr, float);     //! Resolution of the event time of the track computed with the T0AC
DECLARE_SOA_COLUMN(LastTRDCluster, lastTRDCluster, int8_t);  //! Index of the last cluster in the TRD, -1 if no TRD information
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF,                   //! Flag to check if track has a TOF measurement
                           [](float tofSignal) -> bool { return tofSignal > 0; });
// Calibration information

} // namespace tofskims

DECLARE_SOA_TABLE(SkimmedTOFColl, "AOD", "SKIMMEDTOFCOL", //! Table of the skimmed TOF data format. One entry per collision.
                  o2::soa::Index<>,
                  pidtofevtime::EvTimeTOF,
                  pidtofevtime::EvTimeTOFErr,
                  pidtofevtime::EvTimeTOFMult,
                  ft0::T0ACorrected,
                  ft0::T0CCorrected,
                  tofskims::EvTimeT0AC,
                  tofskims::EvTimeT0ACErr,
                  collision::CollisionTime,
                  collision::CollisionTimeRes,
                  pidflags::TOFFlags);

DECLARE_SOA_TABLE(SkimmedTOF, "AOD", "SKIMMEDTOF", //! Table of the skimmed TOF data format. One entry per track.
                  o2::soa::Index<>,
                  tofskims::CollisionId,
                  tofskims::P,
                  tofskims::PtSigned,
                  tofskims::Eta,
                  tofskims::Phi,
                  tofskims::PIDForTracking,
                  track::TOFExpMom,
                  track::Length,
                  track::TOFChi2,
                  track::TPCSignal,
                  pidtofsignal::TOFSignal,
                  pidtofevtime::EvTimeTOF,
                  pidtofevtime::EvTimeTOFErr,
                  tofskims::EvTimeT0AC,
                  tofskims::EvTimeT0ACErr,
                  pidflags::TOFFlags,
                  track::TPCInnerParam,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  track::TPCNClsShared,
                  tofskims::LastTRDCluster,
                  tofskims::HasTOF<pidtofsignal::TOFSignal>,
                  pidflags::IsEvTimeDefined<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOF<pidflags::TOFFlags>,
                  pidflags::IsEvTimeT0AC<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOFT0AC<pidflags::TOFFlags>);

DECLARE_SOA_TABLE(DeltaTOF, "AOD", "DELTATOF", //! Table of the delta TOF data format. One entry per track couple.
                  o2::soa::Index<>,
                  tofskims::CollisionId,
                  tofskims::P,
                  tofskims::DeltaP,
                  tofskims::PtSigned,
                  tofskims::DeltaPt,
                  tofskims::Eta,
                  tofskims::DeltaEta,
                  tofskims::Phi,
                  tofskims::DeltaPhi,
                  tofskims::DeltaTEl,
                  tofskims::DeltaTMu,
                  tofskims::DeltaTPi,
                  tofskims::DeltaTKa,
                  tofskims::DeltaTPr,
                  tofskims::DoubleDelta,
                  tofskims::RefSign,
                  track::Length,
                  track::TOFChi2,
                  track::TPCSignal,
                  pidtofsignal::TOFSignal,
                  pidtofevtime::EvTimeTOF,
                  pidtofevtime::EvTimeTOFErr,
                  pidtofevtime::EvTimeTOFMult,
                  ft0::T0ACorrected,
                  ft0::T0CCorrected,
                  tofskims::EvTimeT0AC,
                  tofskims::EvTimeT0ACErr,
                  collision::CollisionTime,
                  collision::CollisionTimeRes,
                  pidflags::TOFFlags,
                  tofskims::LastTRDCluster);

} // namespace o2::aod

#endif // DPG_TASKS_TOF_TOFSKIMSTABLECREATOR_H_
