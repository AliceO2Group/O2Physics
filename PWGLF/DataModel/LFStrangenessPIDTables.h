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

// Defines TOF PID tables for strangeness.
// Entries calculated per candidate, tables are joinable with v0/cascdata tables.

#ifndef PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
#define PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

namespace o2::aod
{
namespace dautrack
{
// ==== TPC INFORMATION ===
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);     //! track TPC signal
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float); //! Nsigma proton

// ==== TOF INFORMATION ===
DECLARE_SOA_COLUMN(TOFSignal, tofSignal, float); //! track TOF signal
DECLARE_SOA_COLUMN(TOFEvTime, tofEvTime, float); //! track TOF signal
DECLARE_SOA_COLUMN(Length, length, float);       //! track TOF signal
} // namespace dautrack

DECLARE_SOA_TABLE(DauTrackTPCPIDs, "AOD", "DAUTRACKTPCPID", // nsigma table (for analysis)
                  dautrack::TPCSignal, dautrack::TPCNSigmaEl,
                  dautrack::TPCNSigmaPi, dautrack::TPCNSigmaKa,
                  dautrack::TPCNSigmaPr, dautrack::TPCNSigmaHe);
DECLARE_SOA_TABLE(DauTrackTOFPIDs, "AOD", "DAUTRACKTOFPID", // raw table (for posterior TOF calculation)
                  dautrack::TOFSignal, dautrack::TOFEvTime, dautrack::Length);

namespace v0data
{
// ==== TOF INFORMATION ===
// lengths as stored in the AO2D for TOF calculations
DECLARE_SOA_COLUMN(PosTOFLengthToPV, posTOFLengthToPV, float); //! positive track length to PV
DECLARE_SOA_COLUMN(NegTOFLengthToPV, negTOFLengthToPV, float); //! negative track length to PV
DECLARE_SOA_COLUMN(PosTOFSignal, posTOFSignal, float);         //! positive track signal
DECLARE_SOA_COLUMN(NegTOFSignal, negTOFSignal, float);         //! negative track signal
DECLARE_SOA_COLUMN(PosTOFEventTime, posTOFEventTime, float);   //! positive track event time
DECLARE_SOA_COLUMN(NegTOFEventTime, negTOFEventTime, float);   //! negative track event time
DECLARE_SOA_COLUMN(PosTOFLength, posTOFLength, float); //! positive track length, recalculated
DECLARE_SOA_COLUMN(NegTOFLength, negTOFLength, float); //! negative track length, recalculated

// delta-times
DECLARE_SOA_COLUMN(PosTOFDeltaTLaPi, posTOFDeltaTLaPi, float); //! positive track TOFDeltaT from pion <- lambda expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTLaPr, posTOFDeltaTLaPr, float); //! positive track TOFDeltaT from proton <- lambda expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTLaPi, negTOFDeltaTLaPi, float); //! negative track TOFDeltaT from pion <- lambda expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTLaPr, negTOFDeltaTLaPr, float); //! negative track TOFDeltaT from proton <- lambda expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTK0Pi, posTOFDeltaTK0Pi, float); //! positive track TOFDeltaT from pion <- k0short expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTK0Pi, negTOFDeltaTK0Pi, float); //! positive track TOFDeltaT from pion <- k0short expectation

// delta-decay-times (event-time-independent)
DECLARE_SOA_COLUMN(DeltaDecayTimeLambda, deltaDecayTimeLambda, float);         //! delta-decay time estimate from proton/pion from Lambda
DECLARE_SOA_COLUMN(DeltaDecayTimeAntiLambda, deltaDecayTimeAntiLambda, float); //! delta-decay time estimate from pion/proton from AntiLambda
DECLARE_SOA_COLUMN(DeltaDecayTimeK0Short, deltaDecayTimeK0Short, float);       //! delta-decay time estimate from pion/pion from K0Short

// n-sigmas
DECLARE_SOA_COLUMN(TOFNSigmaLaPr, tofNSigmaLaPr, float);           //! positive track NSigma from proton <- lambda expectation
DECLARE_SOA_COLUMN(TOFNSigmaLaPi, tofNSigmaLaPi, float);           //! negative track NSigma from pion <- lambda expectation
DECLARE_SOA_COLUMN(TOFNSigmaALaPr, tofNSigmaALaPr, float);         //! negative track NSigma from proton <- antilambda expectation
DECLARE_SOA_COLUMN(TOFNSigmaALaPi, tofNSigmaALaPi, float);         //! positive track NSigma from pion <- antilambda expectation
DECLARE_SOA_COLUMN(TOFNSigmaK0PiPlus, tofNSigmaK0PiPlus, float);   //! positive track NSigma from pion <- k0short expectation
DECLARE_SOA_COLUMN(TOFNSigmaK0PiMinus, tofNSigmaK0PiMinus, float); //! negative track NSigma from pion <- k0short expectation

// beta values
DECLARE_SOA_COLUMN(TofBetaLambda, tofBetaLambda, float);         //! beta value with Lambda hypothesis
DECLARE_SOA_COLUMN(TofBetaAntiLambda, tofBetaAntiLambda, float); //! beta value with AntiLambda hypothesis
DECLARE_SOA_COLUMN(TofBetaK0Short, tofBetaK0Short, float);       //! beta value with K0Short hypothesis

// debug quantities
DECLARE_SOA_COLUMN(V0LifetimeLambda, v0LifetimeLambda, float);   //! lifetime of V0 assuming lambda mass (ps)
DECLARE_SOA_COLUMN(V0LifetimeK0Short, v0LifetimeK0Short, float); //! lifetime of V0 assuming K0 mass (ps)
DECLARE_SOA_COLUMN(PosLifetimePr, posLifetimePr, float);         //! lifetime (to TOF) of pos prong assuming proton mass (ps)
DECLARE_SOA_COLUMN(PosLifetimePi, posLifetimePi, float);         //! lifetime (to TOF) of pos prong assuming pion mass (ps)
DECLARE_SOA_COLUMN(NegLifetimePr, negLifetimePr, float);         //! lifetime (to TOF) of neg prong assuming proton mass (ps)
DECLARE_SOA_COLUMN(NegLifetimePi, negLifetimePi, float);         //! lifetime (to TOF) of neg prong assuming pion mass (ps)
} // namespace v0data

// /-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/
// DEPRECATED - DO NOT USE - KEPT FOR BACKWARDS COMPATIBILITY, TO BE REMOVED
DECLARE_SOA_TABLE(V0TOFs, "AOD", "V0TOF", // raw information table (for debug, etc)
                  v0data::PosTOFLengthToPV, v0data::NegTOFLengthToPV,
                  v0data::PosTOFSignal, v0data::NegTOFSignal,
                  v0data::PosTOFEventTime, v0data::NegTOFEventTime);
// /-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/

DECLARE_SOA_TABLE(V0TOFPIDs, "AOD", "V0TOFPID", // processed info table (for analysis)
                  v0data::PosTOFDeltaTLaPi, v0data::PosTOFDeltaTLaPr,
                  v0data::NegTOFDeltaTLaPi, v0data::NegTOFDeltaTLaPr,
                  v0data::PosTOFDeltaTK0Pi, v0data::NegTOFDeltaTK0Pi,
                  v0data::DeltaDecayTimeLambda,
                  v0data::DeltaDecayTimeAntiLambda,
                  v0data::DeltaDecayTimeK0Short);

DECLARE_SOA_TABLE(V0TOFDebugs, "AOD", "V0TOFDEBUG", // table with intermediate information solely for debugging
                  v0data::V0LifetimeLambda, v0data::V0LifetimeK0Short,
                  v0data::PosLifetimePr, v0data::PosLifetimePi,
                  v0data::NegLifetimePr, v0data::NegLifetimePi);

DECLARE_SOA_TABLE(V0TOFBetas, "AOD", "V0TOFBETA", // processed info table (for analysis)
                  v0data::TofBetaLambda, v0data::TofBetaAntiLambda, v0data::TofBetaK0Short);

DECLARE_SOA_TABLE(V0TOFNSigmas, "AOD", "V0TOFNSIGMA", // processed NSigma table (for analysis)
                  v0data::TOFNSigmaLaPr, v0data::TOFNSigmaLaPi,
                  v0data::TOFNSigmaALaPr, v0data::TOFNSigmaALaPi,
                  v0data::TOFNSigmaK0PiPlus, v0data::TOFNSigmaK0PiMinus);

namespace cascdata
{
// lengths as stored in the AO2D for TOF calculations
DECLARE_SOA_COLUMN(PosTOFLengthToPV, posTOFLengthToPV, float);   //! positive track length
DECLARE_SOA_COLUMN(NegTOFLengthToPV, negTOFLengthToPV, float);   //! negative track length
DECLARE_SOA_COLUMN(BachTOFLengthToPV, bachTOFLengthToPV, float); //! bachelor track length
DECLARE_SOA_COLUMN(PosTOFSignal, posTOFSignal, float);           //! positive track signal
DECLARE_SOA_COLUMN(NegTOFSignal, negTOFSignal, float);           //! negative track signal
DECLARE_SOA_COLUMN(BachTOFSignal, bachTOFSignal, float);         //! bachelor track signal
DECLARE_SOA_COLUMN(PosTOFEventTime, posTOFEventTime, float);     //! positive track event time
DECLARE_SOA_COLUMN(NegTOFEventTime, negTOFEventTime, float);     //! negative track event time
DECLARE_SOA_COLUMN(BachTOFEventTime, bachTOFEventTime, float);   //! bachelor track event time

// delta-times
DECLARE_SOA_COLUMN(PosTOFDeltaTXiPi, posTOFDeltaTXiPi, float);   //! positive track TOFDeltaT from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTXiPr, posTOFDeltaTXiPr, float);   //! positive track TOFDeltaT from proton <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTXiPi, negTOFDeltaTXiPi, float);   //! negative track TOFDeltaT from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTXiPr, negTOFDeltaTXiPr, float);   //! negative track TOFDeltaT from proton <- lambda <- xi expectation
DECLARE_SOA_COLUMN(BachTOFDeltaTXiPi, bachTOFDeltaTXiPi, float); //! bachelor track TOFDeltaT from pion <- xi expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTOmPi, posTOFDeltaTOmPi, float);   //! positive track TOFDeltaT from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTOmPr, posTOFDeltaTOmPr, float);   //! positive track TOFDeltaT from proton <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTOmPi, negTOFDeltaTOmPi, float);   //! negative track TOFDeltaT from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTOmPr, negTOFDeltaTOmPr, float);   //! negative track TOFDeltaT from proton <- lambda <- omega expectation
DECLARE_SOA_COLUMN(BachTOFDeltaTOmKa, bachTOFDeltaTOmKa, float); //! bachelor track TOFDeltaT from kaon <- omega expectation

// n-sigmas
DECLARE_SOA_COLUMN(TOFNSigmaXiLaPi, tofNSigmaXiLaPi, float); //! meson track NSigma from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(TOFNSigmaXiLaPr, tofNSigmaXiLaPr, float); //! baryon track NSigma from proton <- lambda <- xi expectation
DECLARE_SOA_COLUMN(TOFNSigmaXiPi, tofNSigmaXiPi, float);     //! bachelor track NSigma from pion <- xi expectation
DECLARE_SOA_COLUMN(TOFNSigmaOmLaPi, tofNSigmaOmLaPi, float); //! meson track NSigma from pion <- lambda <- om expectation
DECLARE_SOA_COLUMN(TOFNSigmaOmLaPr, tofNSigmaOmLaPr, float); //! baryon track NSigma from proton <- lambda <- om expectation
DECLARE_SOA_COLUMN(TOFNSigmaOmKa, tofNSigmaOmKa, float);     //! bachelor track NSigma from kaon <- om expectation
} // namespace cascdata

// /-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/
// DEPRECATED - DO NOT USE - KEPT FOR BACKWARDS COMPATIBILITY, TO BE REMOVED
DECLARE_SOA_TABLE(CascTOFs, "AOD", "CascTOF", // raw information table (for debug, etc)
                  cascdata::PosTOFLengthToPV, cascdata::NegTOFLengthToPV, cascdata::BachTOFLengthToPV,
                  cascdata::PosTOFSignal, cascdata::NegTOFSignal, cascdata::BachTOFSignal,
                  cascdata::PosTOFEventTime, cascdata::NegTOFEventTime, cascdata::BachTOFEventTime);
// /-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/-|-\-|-/

DECLARE_SOA_TABLE(CascTOFPIDs, "AOD", "CASCTOFPID", // processed information for analysis
                  cascdata::PosTOFDeltaTXiPi, cascdata::PosTOFDeltaTXiPr,
                  cascdata::NegTOFDeltaTXiPi, cascdata::NegTOFDeltaTXiPr,
                  cascdata::BachTOFDeltaTXiPi,
                  cascdata::PosTOFDeltaTOmPi, cascdata::PosTOFDeltaTOmPr,
                  cascdata::NegTOFDeltaTOmPi, cascdata::NegTOFDeltaTOmPr,
                  cascdata::BachTOFDeltaTOmKa);
DECLARE_SOA_TABLE(CascTOFNSigmas, "AOD", "CascTOFNSigmas", // Nsigmas for cascades
                  cascdata::TOFNSigmaXiLaPi, cascdata::TOFNSigmaXiLaPr, cascdata::TOFNSigmaXiPi,
                  cascdata::TOFNSigmaOmLaPi, cascdata::TOFNSigmaOmLaPr, cascdata::TOFNSigmaOmKa);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
