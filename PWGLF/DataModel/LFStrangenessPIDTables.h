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

//**********************************************************************
// Defines TOF PID tables for strangeness.
// Entries calculated per candidate, tables are joinable with v0/cascdata tables.
//**********************************************************************

//**********************************************************************
// Nota bene: when using, do not check track.hasTOF! That conditional may not match
// the calculation of strangeness TOF, which requires e.g. a successful calculation
// of the collision time for the reassociated collision
//**********************************************************************

#ifndef PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
#define PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace dautrack
{
// ==== define packing helpers ===
namespace packing
{
// define variables for packing
static constexpr int nbins = (1 << 8 * sizeof(int8_t)) - 2;
static constexpr int8_t overflowBin = nbins >> 1;
static constexpr int8_t underflowBin = -(nbins >> 1);
static constexpr float binned_max = 6.35;
static constexpr float binned_min = -6.35;
static constexpr float bin_width = (binned_max - binned_min) / nbins;
static constexpr float underflow_return = -100.0f;
static constexpr float overflow_return = +100.0f;

// define helper function to do packing
int8_t packInInt8(float nSigma)
{
  // calculate
  if (nSigma <= binned_min)
    return underflowBin;
  if (nSigma >= binned_max)
    return overflowBin;
  if (nSigma >= 0) {
    return static_cast<int8_t>((nSigma / bin_width) + 0.5f);
  }
  // automatic: this is the case in which nSigma < 0
  return static_cast<int8_t>((nSigma / bin_width) - 0.5f);
}

// define helper function to do unpacking
float unpackInt8(int8_t nSigma)
{
  if (nSigma == underflowBin) {
    return underflow_return;
  }
  if (nSigma == overflowBin) {
    return overflow_return;
  }
  return bin_width * nSigma;
}

} // namespace packing
} // namespace dautrack

namespace dautrack_legacy
{
// ==== LEGACY TPC INFORMATION (full size tables) ===
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float); //! Nsigma proton
} // namespace dautrack_legacy

namespace dautrack
{
// ==== COMPACT TPC INFORMATION (full size tables) ===
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                  //! track TPC signal
DECLARE_SOA_COLUMN(PackedTPCNSigmaEl, packedTpcNSigmaEl, int8_t); //! Nsigma proton
DECLARE_SOA_COLUMN(PackedTPCNSigmaPi, packedTpcNSigmaPi, int8_t); //! Nsigma proton
DECLARE_SOA_COLUMN(PackedTPCNSigmaKa, packedTpcNSigmaKa, int8_t); //! Nsigma proton
DECLARE_SOA_COLUMN(PackedTPCNSigmaPr, packedTpcNSigmaPr, int8_t); //! Nsigma proton

DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaEl, tpcNSigmaEl, //! unpacked TPC nsigma
                           [](int8_t nsigma_packed) -> float { return o2::aod::dautrack::packing::unpackInt8(nsigma_packed); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi, tpcNSigmaPi, //! unpacked TPC nsigma
                           [](int8_t nsigma_packed) -> float { return o2::aod::dautrack::packing::unpackInt8(nsigma_packed); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa, tpcNSigmaKa, //! unpacked TPC nsigma
                           [](int8_t nsigma_packed) -> float { return o2::aod::dautrack::packing::unpackInt8(nsigma_packed); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr, //! unpacked TPC nsigma
                           [](int8_t nsigma_packed) -> float { return o2::aod::dautrack::packing::unpackInt8(nsigma_packed); });

// ==== TOF INFORMATION ===
DECLARE_SOA_INDEX_COLUMN(DauTrackExtra, dauTrackExtra); //! point to daughter this TOF info belongs to
DECLARE_SOA_INDEX_COLUMN(StraCollision, straCollision); //! point to collision associated with this track (not the V0/Casc)
DECLARE_SOA_COLUMN(TOFSignal, tofSignal, float);        //! track TOF signal
DECLARE_SOA_COLUMN(TOFEvTime, tofEvTime, float);        //! event time
DECLARE_SOA_COLUMN(TOFEvTimeErr, tofEvTimeErr, float);  //! event time error for TOF
DECLARE_SOA_COLUMN(Length, length, float);              //! track length (to assigned PV)
DECLARE_SOA_COLUMN(TOFExpMom, tofExpMom, float);        //! tof Exp Mom (to assigned PV)

// dynamics with expected times
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpTimeEl, tofExpTimeEl, //! Expected time for the track to reach the TOF under the electron hypothesis
                           [](float length, float tofExpMom) -> float {
                             constexpr float massSquared = o2::constants::physics::MassElectron * o2::constants::physics::MassElectron;
                             return o2::framework::pid::tof::MassToExpTime(tofExpMom, length, massSquared);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(TOFExpTimePi, tofExpTimePi, //! Expected time for the track to reach the TOF under the pion hypothesis
                           [](float length, float tofExpMom) -> float {
                             constexpr float massSquared = o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged;
                             return o2::framework::pid::tof::MassToExpTime(tofExpMom, length, massSquared);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(TOFExpTimeKa, tofExpTimeKa, //! Expected time for the track to reach the TOF under the kaon hypothesis
                           [](float length, float tofExpMom) -> float {
                             constexpr float massSquared = o2::constants::physics::MassKaonCharged * o2::constants::physics::MassKaonCharged;
                             return o2::framework::pid::tof::MassToExpTime(tofExpMom, length, massSquared);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(TOFExpTimePr, tofExpTimePr, //! Expected time for the track to reach the TOF under the proton hypothesis
                           [](float length, float tofExpMom) -> float {
                             constexpr float massSquared = o2::constants::physics::MassProton * o2::constants::physics::MassProton;
                             return o2::framework::pid::tof::MassToExpTime(tofExpMom, length, massSquared);
                           });

} // namespace dautrack

DECLARE_SOA_TABLE(DauTrackTPCPIDs_000, "AOD", "DAUTRACKTPCPID", // nsigma table (for analysis)
                  dautrack::TPCSignal, dautrack_legacy::TPCNSigmaEl,
                  dautrack_legacy::TPCNSigmaPi, dautrack_legacy::TPCNSigmaKa,
                  dautrack_legacy::TPCNSigmaPr, dautrack_legacy::TPCNSigmaHe);

DECLARE_SOA_TABLE_VERSIONED(DauTrackTPCPIDs_001, "AOD", "DAUTRACKTPCPID", 1, // nsigma table (for analysis)
                            dautrack::TPCSignal,
                            dautrack::PackedTPCNSigmaEl, dautrack::PackedTPCNSigmaPi,
                            dautrack::PackedTPCNSigmaKa, dautrack::PackedTPCNSigmaPr,
                            dautrack::TPCNSigmaEl<dautrack::PackedTPCNSigmaEl>,
                            dautrack::TPCNSigmaPi<dautrack::PackedTPCNSigmaPi>,
                            dautrack::TPCNSigmaKa<dautrack::PackedTPCNSigmaKa>,
                            dautrack::TPCNSigmaPr<dautrack::PackedTPCNSigmaPr>);

using DauTrackTPCPIDs = DauTrackTPCPIDs_001; // second gen: packed Nsigma, no He

DECLARE_SOA_TABLE(DauTrackTOFPIDs_000, "AOD", "DAUTRACKTOFPID", // raw table (for posterior TOF calculation)
                  dautrack::TOFSignal, dautrack::TOFEvTime, dautrack::Length);

DECLARE_SOA_TABLE_VERSIONED(DauTrackTOFPIDs_001, "AOD", "DAUTRACKTOFPID", 1, // raw table (for posterior TOF calculation)
                            o2::soa::Index<>,
                            dautrack::StraCollisionId, dautrack::DauTrackExtraId,
                            dautrack::TOFSignal, dautrack::TOFEvTime,
                            dautrack::Length, dautrack::TOFExpMom,
                            dautrack::TOFExpTimeEl<dautrack::Length, dautrack::TOFExpMom>,
                            dautrack::TOFExpTimePi<dautrack::Length, dautrack::TOFExpMom>,
                            dautrack::TOFExpTimeKa<dautrack::Length, dautrack::TOFExpMom>,
                            dautrack::TOFExpTimePr<dautrack::Length, dautrack::TOFExpMom>);

DECLARE_SOA_TABLE_VERSIONED(DauTrackTOFPIDs_002, "AOD", "DAUTRACKTOFPID", 2, // raw table (for posterior TOF calculation)
                            o2::soa::Index<>,
                            dautrack::StraCollisionId, dautrack::DauTrackExtraId,
                            dautrack::TOFSignal, dautrack::TOFEvTime, dautrack::TOFEvTimeErr,
                            dautrack::Length, dautrack::TOFExpMom,
                            dautrack::TOFExpTimeEl<dautrack::Length, dautrack::TOFExpMom>,
                            dautrack::TOFExpTimePi<dautrack::Length, dautrack::TOFExpMom>,
                            dautrack::TOFExpTimeKa<dautrack::Length, dautrack::TOFExpMom>,
                            dautrack::TOFExpTimePr<dautrack::Length, dautrack::TOFExpMom>);

using DauTrackTOFPIDs = DauTrackTOFPIDs_002; // second gen: with collision Id, with TOFExpMom

namespace v0data
{
// define constants for NSigma operation
constexpr float kNoTOFValue = -1e+6;
const float kEpsilon = 1e-4;

// ==== TOF INFORMATION ===
// lengths as stored in the AO2D for TOF calculations
DECLARE_SOA_COLUMN(PosTOFLengthToPV, posTOFLengthToPV, float); //! positive track length to PV
DECLARE_SOA_COLUMN(NegTOFLengthToPV, negTOFLengthToPV, float); //! negative track length to PV
DECLARE_SOA_COLUMN(PosTOFSignal, posTOFSignal, float);         //! positive track signal
DECLARE_SOA_COLUMN(NegTOFSignal, negTOFSignal, float);         //! negative track signal
DECLARE_SOA_COLUMN(PosTOFEventTime, posTOFEventTime, float);   //! positive track event time
DECLARE_SOA_COLUMN(NegTOFEventTime, negTOFEventTime, float);   //! negative track event time
DECLARE_SOA_COLUMN(PosTOFLength, posTOFLength, float);         //! positive track length, recalculated
DECLARE_SOA_COLUMN(NegTOFLength, negTOFLength, float);         //! negative track length, recalculated

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

// for wrong hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElPosFromPhoton, tofNSigmaElPosFromPhoton, float); //! n sigma of positive track from photon conversion under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElNegFromPhoton, tofNSigmaElNegFromPhoton, float); //! n sigma of negative track from photon conversion under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiPosFromPhoton, tofNSigmaPiPosFromPhoton, float); //! n sigma of positive track from photon conversion under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiNegFromPhoton, tofNSigmaPiNegFromPhoton, float); //! n sigma of negative track from photon conversion under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaPosFromPhoton, tofNSigmaKaPosFromPhoton, float); //! n sigma of positive track from photon conversion under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaNegFromPhoton, tofNSigmaKaNegFromPhoton, float); //! n sigma of negative track from photon conversion under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrPosFromPhoton, tofNSigmaPrPosFromPhoton, float); //! n sigma of positive track from photon conversion under proton hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrNegFromPhoton, tofNSigmaPrNegFromPhoton, float); //! n sigma of negative track from photon conversion under proton hypothesis

DECLARE_SOA_COLUMN(TOFNSigmaElPosFromK0S, tofNSigmaElPosFromK0S, float); //! n sigma of positive track from K0S under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElNegFromK0S, tofNSigmaElNegFromK0S, float); //! n sigma of negative track from K0S under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiPosFromK0S, tofNSigmaPiPosFromK0S, float); //! n sigma of positive track from K0S under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiNegFromK0S, tofNSigmaPiNegFromK0S, float); //! n sigma of negative track from K0S under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaPosFromK0S, tofNSigmaKaPosFromK0S, float); //! n sigma of positive track from K0S under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaNegFromK0S, tofNSigmaKaNegFromK0S, float); //! n sigma of negative track from K0S under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrPosFromK0S, tofNSigmaPrPosFromK0S, float); //! n sigma of positive track from K0S under proton hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrNegFromK0S, tofNSigmaPrNegFromK0S, float); //! n sigma of negative track from K0S under proton hypothesis

DECLARE_SOA_COLUMN(TOFNSigmaElPosFromLambda, tofNSigmaElPosFromLambda, float); //! n sigma of positive track from Lambda under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElNegFromLambda, tofNSigmaElNegFromLambda, float); //! n sigma of negative track from Lambda under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiPosFromLambda, tofNSigmaPiPosFromLambda, float); //! n sigma of positive track from Lambda under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiNegFromLambda, tofNSigmaPiNegFromLambda, float); //! n sigma of negative track from Lambda under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaPosFromLambda, tofNSigmaKaPosFromLambda, float); //! n sigma of positive track from Lambda under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaNegFromLambda, tofNSigmaKaNegFromLambda, float); //! n sigma of negative track from Lambda under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrPosFromLambda, tofNSigmaPrPosFromLambda, float); //! n sigma of positive track from Lambda under proton hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrNegFromLambda, tofNSigmaPrNegFromLambda, float); //! n sigma of negative track from Lambda under proton hypothesis

// dynamics to replace hasTOF (note: that condition does not match track hasTOF!)
// note: only single hypothesis check necessary; other hypotheses will always be valid
DECLARE_SOA_DYNAMIC_COLUMN(PositiveHasTOF, positiveHasTOF, //! positive daughter TOF calculation valid
                           [](float TOFNSigmaLaPr) -> bool {
                             bool returnStatus = true;
                             if (std::abs(TOFNSigmaLaPr - kNoTOFValue) < kEpsilon) {
                               returnStatus = false;
                             }
                             return returnStatus;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeHasTOF, negativeHasTOF, //! negative daughter TOF calculation valid
                           [](float TOFNSigmaALaPr) -> bool {
                             bool returnStatus = true;
                             if (std::abs(TOFNSigmaALaPr - kNoTOFValue) < kEpsilon) {
                               returnStatus = false;
                             }
                             return returnStatus;
                           });

// dynamics based on n-sigmas with use-only-if-tof-present logic
DECLARE_SOA_DYNAMIC_COLUMN(TofLambdaCompatibility, tofLambdaCompatibility, //! compatibility with being lambda, checked only if TOF present. Argument: number of sigmas
                           [](float tofNSigmaLaPr, float tofNSigmaLaPi, float nsigma) -> float {
                             bool compatible = true;
                             if (std::abs(tofNSigmaLaPr - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaLaPr) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             if (std::abs(tofNSigmaLaPi - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaLaPi) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             return compatible;
                           });

// dynamics based on n-sigmas with use-only-if-tof-present logic
DECLARE_SOA_DYNAMIC_COLUMN(TofAntiLambdaCompatibility, tofAntiLambdaCompatibility, //! compatibility with being lambda, checked only if TOF present. Argument: number of sigmas
                           [](float tofNSigmaALaPr, float tofNSigmaALaPi, float nsigma) -> float {
                             bool compatible = true;
                             if (std::abs(tofNSigmaALaPr - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaALaPr) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             if (std::abs(tofNSigmaALaPi - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaALaPi) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             return compatible;
                           });

// dynamics based on n-sigmas with use-only-if-tof-present logic
DECLARE_SOA_DYNAMIC_COLUMN(TofK0ShortCompatibility, tofK0ShortCompatibility, //! compatibility with being lambda, checked only if TOF present. Argument: number of sigmas
                           [](float tofNSigmaK0PiPlus, float tofNSigmaK0PiMinus, float nsigma) -> float {
                             bool compatible = true;
                             if (std::abs(tofNSigmaK0PiPlus - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaK0PiPlus) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             if (std::abs(tofNSigmaK0PiMinus - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaK0PiMinus) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             return compatible;
                           });

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
                  v0data::TOFNSigmaK0PiPlus, v0data::TOFNSigmaK0PiMinus,
                  v0data::PositiveHasTOF<v0data::TOFNSigmaLaPr>,
                  v0data::NegativeHasTOF<v0data::TOFNSigmaALaPr>,
                  v0data::TofLambdaCompatibility<v0data::TOFNSigmaLaPr, v0data::TOFNSigmaLaPi>,
                  v0data::TofAntiLambdaCompatibility<v0data::TOFNSigmaALaPr, v0data::TOFNSigmaALaPi>,
                  v0data::TofK0ShortCompatibility<v0data::TOFNSigmaK0PiPlus, v0data::TOFNSigmaK0PiMinus>);

DECLARE_SOA_TABLE(V0TOFNSigmasAll, "AOD", "V0TOFNSIGMAALL", // processed NSigma table (for analysis) including wrong hypothesis
                  v0data::TOFNSigmaElPosFromPhoton, v0data::TOFNSigmaElPosFromK0S, v0data::TOFNSigmaElPosFromLambda,
                  v0data::TOFNSigmaElNegFromPhoton, v0data::TOFNSigmaElNegFromK0S, v0data::TOFNSigmaElNegFromLambda,
                  v0data::TOFNSigmaPiPosFromPhoton, v0data::TOFNSigmaPiPosFromK0S, v0data::TOFNSigmaPiPosFromLambda,
                  v0data::TOFNSigmaPiNegFromPhoton, v0data::TOFNSigmaPiNegFromK0S, v0data::TOFNSigmaPiNegFromLambda,
                  v0data::TOFNSigmaKaPosFromPhoton, v0data::TOFNSigmaKaPosFromK0S, v0data::TOFNSigmaKaPosFromLambda,
                  v0data::TOFNSigmaKaNegFromPhoton, v0data::TOFNSigmaKaNegFromK0S, v0data::TOFNSigmaKaNegFromLambda,
                  v0data::TOFNSigmaPrPosFromPhoton, v0data::TOFNSigmaPrPosFromK0S, v0data::TOFNSigmaPrPosFromLambda,
                  v0data::TOFNSigmaPrNegFromPhoton, v0data::TOFNSigmaPrNegFromK0S, v0data::TOFNSigmaPrNegFromLambda);

namespace cascdata
{
// define constants for NSigma operation
const float kNoTOFValue = -1e+6;
const float kEpsilon = 1e-4;

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

// for wrong hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElFromLambdaFromXi, tofNSigmaElFromLambdaFromXi, float);       //! nigma of positive track from Lambda from Xi under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElFromXi, tofNSigmaElFromXi, float);                           //! nigma of bachelor track from Xi under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElFromLambdaFromOmega, tofNSigmaElFromLambdaFromOmega, float); //! nigma of positive track from Lambda from Omega under electron hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaElFromOmega, tofNSigmaElFromOmega, float);                     //! nigma of bachelor track from Omega under electron hypothesis

DECLARE_SOA_COLUMN(TOFNSigmaPiFromLambdaFromXi, tofNSigmaPiFromLambdaFromXi, float);       //! nigma of positive track from Lambda from Xi under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiFromXi, tofNSigmaPiFromXi, float);                           //! nigma of bachelor track from Xi under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiFromLambdaFromOmega, tofNSigmaPiFromLambdaFromOmega, float); //! nigma of positive track from Lambda from Omega under pion hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPiFromOmega, tofNSigmaPiFromOmega, float);                     //! nigma of bachelor track from Omega under pion hypothesis

DECLARE_SOA_COLUMN(TOFNSigmaKaFromLambdaFromXi, tofNSigmaKaFromLambdaFromXi, float);       //! nigma of positive track from Lambda from Xi under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaFromXi, tofNSigmaKaFromXi, float);                           //! nigma of bachelor track from Xi under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaFromLambdaFromOmega, tofNSigmaKaFromLambdaFromOmega, float); //! nigma of positive track from Lambda from Omega under kaon hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaKaFromOmega, tofNSigmaKaFromOmega, float);                     //! nigma of bachelor track from Omega under kaon hypothesis

DECLARE_SOA_COLUMN(TOFNSigmaPrFromLambdaFromXi, tofNSigmaPrFromLambdaFromXi, float);       //! nigma of positive track from Lambda from Xi under proton hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrFromXi, tofNSigmaPrFromXi, float);                           //! nigma of bachelor track from Xi under proton hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrFromLambdaFromOmega, tofNSigmaPrFromLambdaFromOmega, float); //! nigma of positive track from Lambda from Omega under proton hypothesis
DECLARE_SOA_COLUMN(TOFNSigmaPrFromOmega, tofNSigmaPrFromOmega, float);                     //! nigma of bachelor track from Omega under proton hypothesis

// dynamics to replace hasTOF (note: that condition does not match track hasTOF!)
// note: only single hypothesis check necessary; other hypotheses will always be valid
DECLARE_SOA_DYNAMIC_COLUMN(PositiveHasTOF, positiveHasTOF, //! positive daughter TOF calculation valid
                           [](float PosTOFDeltaTXiPr) -> bool {
                             bool returnStatus = true;
                             if (std::abs(PosTOFDeltaTXiPr - kNoTOFValue) < kEpsilon) {
                               returnStatus = false;
                             }
                             return returnStatus;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeHasTOF, negativeHasTOF, //! positive daughter TOF calculation valid
                           [](float NegTOFDeltaTXiPr) -> bool {
                             bool returnStatus = true;
                             if (std::abs(NegTOFDeltaTXiPr - kNoTOFValue) < kEpsilon) {
                               returnStatus = false;
                             }
                             return returnStatus;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(BachelorHasTOF, bachelorHasTOF, //! bachelor daughter TOF calculation valid
                           [](float BachTOFDeltaTXiPi) -> bool {
                             bool returnStatus = true;
                             if (std::abs(BachTOFDeltaTXiPi - kNoTOFValue) < kEpsilon) {
                               returnStatus = false;
                             }
                             return returnStatus;
                           });

// dynamics based on n-sigmas with use-only-if-tof-present logic
DECLARE_SOA_DYNAMIC_COLUMN(TofXiCompatibility, tofXiCompatibility, //! compatibility with being lambda, checked only if TOF present. Argument: number of sigmas
                           [](float tofNSigmaXiLaPr, float tofNSigmaXiLaPi, float tofNSigmaXiPi, float nsigma) -> float {
                             bool compatible = true;
                             if (std::abs(tofNSigmaXiLaPr - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaXiLaPr) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             if (std::abs(tofNSigmaXiLaPi - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaXiLaPi) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             if (std::abs(tofNSigmaXiPi - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaXiPi) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             return compatible;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(TofOmegaCompatibility, tofOmegaCompatibility, //! compatibility with being lambda, checked only if TOF present. Argument: number of sigmas
                           [](float tofNSigmaOmLaPr, float tofNSigmaOmLaPi, float tofNSigmaOmKa, float nsigma) -> float {
                             bool compatible = true;
                             if (std::abs(tofNSigmaOmLaPr - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaOmLaPr) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             if (std::abs(tofNSigmaOmLaPi - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaOmLaPi) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             if (std::abs(tofNSigmaOmKa - kNoTOFValue) > kEpsilon && std::abs(tofNSigmaOmKa) > nsigma) {
                               compatible = false; // reject only if info present and incompatible
                             }
                             return compatible;
                           });

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
                  cascdata::TOFNSigmaOmLaPi, cascdata::TOFNSigmaOmLaPr, cascdata::TOFNSigmaOmKa,
                  cascdata::PositiveHasTOF<cascdata::PosTOFDeltaTXiPr>,
                  cascdata::NegativeHasTOF<cascdata::NegTOFDeltaTXiPr>,
                  cascdata::BachelorHasTOF<cascdata::BachTOFDeltaTXiPi>,
                  cascdata::TofXiCompatibility<cascdata::TOFNSigmaXiLaPr, cascdata::TOFNSigmaXiLaPi, cascdata::TOFNSigmaXiPi>,
                  cascdata::TofOmegaCompatibility<cascdata::TOFNSigmaOmLaPr, cascdata::TOFNSigmaOmLaPi, cascdata::TOFNSigmaOmKa>);

DECLARE_SOA_TABLE(CascTOFNSigmasAll, "AOD", "CascTOFNSigmasAll", // Nsigmas for cascades including wrong hypothesis
                  cascdata::TOFNSigmaElFromLambdaFromXi, cascdata::TOFNSigmaElFromXi, cascdata::TOFNSigmaElFromLambdaFromOmega, cascdata::TOFNSigmaElFromOmega,
                  cascdata::TOFNSigmaPiFromLambdaFromXi, cascdata::TOFNSigmaPiFromXi, cascdata::TOFNSigmaPiFromLambdaFromOmega, cascdata::TOFNSigmaPiFromOmega,
                  cascdata::TOFNSigmaKaFromLambdaFromXi, cascdata::TOFNSigmaKaFromXi, cascdata::TOFNSigmaKaFromLambdaFromOmega, cascdata::TOFNSigmaKaFromOmega,
                  cascdata::TOFNSigmaPrFromLambdaFromXi, cascdata::TOFNSigmaPrFromXi, cascdata::TOFNSigmaPrFromLambdaFromOmega, cascdata::TOFNSigmaPrFromOmega);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
