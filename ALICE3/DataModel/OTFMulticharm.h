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
/// \file   OTFMulticharm.h
/// \author David Dobrigkeit Chinellato
/// \author Jesper Karlsson Gumprecht
/// \since  05/08/2024
/// \brief  Set of tables for the ALICE3 multi-charm information
///

#ifndef ALICE3_DATAMODEL_OTFMULTICHARM_H_
#define ALICE3_DATAMODEL_OTFMULTICHARM_H_

// O2 includes
#include "ALICE3/DataModel/OTFStrangeness.h"

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace otfmulticharm
{
DECLARE_SOA_INDEX_COLUMN_FULL(Cascade, cascade, int, UpgradeCascades, "_Cascade");
DECLARE_SOA_INDEX_COLUMN_FULL(XiCPion1, xiCPion1, int, Tracks, "_Pi1XiC");
DECLARE_SOA_INDEX_COLUMN_FULL(XiCPion2, xiCPion2, int, Tracks, "_Pi2XiC");
DECLARE_SOA_INDEX_COLUMN_FULL(XiCCPion, xiCCPion, int, Tracks, "_PiXiCC");

DECLARE_SOA_COLUMN(XicMass, xicMass, float);
DECLARE_SOA_COLUMN(XiccMass, xiccMass, float);
DECLARE_SOA_COLUMN(LUTConfigId, lutConfigId, int); //! Index for LUT configuration

// kine vars
DECLARE_SOA_COLUMN(XiccPt, xiccPt, float);
DECLARE_SOA_COLUMN(XiccEta, xiccEta, float);
DECLARE_SOA_COLUMN(XicPt, xicPt, float);
DECLARE_SOA_COLUMN(XicEta, xicEta, float);

// topo vars
DECLARE_SOA_COLUMN(XiDCAz, xiDCAz, float);
DECLARE_SOA_COLUMN(XiDCAxy, xiDCAxy, float);
DECLARE_SOA_COLUMN(XicDauDCA, xicDauDCA, float);
DECLARE_SOA_COLUMN(XicDCAxy, xicDCAxy, float);
DECLARE_SOA_COLUMN(XicDCAz, xicDCAz, float);
DECLARE_SOA_COLUMN(XiccDauDCA, xiccDauDCA, float);
DECLARE_SOA_COLUMN(XiccDCAxy, xiccDCAxy, float);
DECLARE_SOA_COLUMN(XiccDCAz, xiccDCAz, float);

DECLARE_SOA_COLUMN(BachDCAxy, bachDCAxy, float);
DECLARE_SOA_COLUMN(BachDCAz, bachDCAz, float);
DECLARE_SOA_COLUMN(PosDCAxy, posDCAxy, float);
DECLARE_SOA_COLUMN(PosDCAz, posDCAz, float);
DECLARE_SOA_COLUMN(NegDCAxy, negDCAxy, float);
DECLARE_SOA_COLUMN(NegDCAz, negDCAz, float);

DECLARE_SOA_COLUMN(Pi1cDCAxy, pi1cDCAxy, float);
DECLARE_SOA_COLUMN(Pi1cDCAz, pi1cDCAz, float);
DECLARE_SOA_COLUMN(Pi2cDCAxy, pi2cDCAxy, float);
DECLARE_SOA_COLUMN(Pi2cDCAz, pi2cDCAz, float);
DECLARE_SOA_COLUMN(PiccDCAxy, piccDCAxy, float);
DECLARE_SOA_COLUMN(PiccDCAz, piccDCAz, float);

// Lengths
DECLARE_SOA_COLUMN(XicDecayRadius2D, xicDecayRadius2D, float);
DECLARE_SOA_COLUMN(XiccDecayRadius2D, xiccDecayRadius2D, float);
DECLARE_SOA_COLUMN(XicProperLength, xicProperLength, float);
DECLARE_SOA_COLUMN(XicDistanceFromPV, xicDistanceFromPV, float);
DECLARE_SOA_COLUMN(XiccProperLength, xiccProperLength, float);

// PID
DECLARE_SOA_COLUMN(Pi1cTofDeltaInner, pi1cTofDeltaInner, float);
DECLARE_SOA_COLUMN(Pi1cTofNSigmaInner, pi1cTofNSigmaInner, float);
DECLARE_SOA_COLUMN(Pi1cTofDeltaOuter, pi1cTofDeltaOuter, float);
DECLARE_SOA_COLUMN(Pi1cTofNSigmaOuter, pi1cTofNSigmaOuter, float);
DECLARE_SOA_COLUMN(Pi1cHasRichSignal, pi1cHasRichSignal, bool);
DECLARE_SOA_COLUMN(Pi1cRichNSigma, pi1cRichNSigma, float);
DECLARE_SOA_COLUMN(Pi1cPdgCode, pi1cPdgCode, int);

DECLARE_SOA_COLUMN(Pi2cTofDeltaInner, pi2cTofDeltaInner, float);
DECLARE_SOA_COLUMN(Pi2cTofNSigmaInner, pi2cTofNSigmaInner, float);
DECLARE_SOA_COLUMN(Pi2cTofDeltaOuter, pi2cTofDeltaOuter, float);
DECLARE_SOA_COLUMN(Pi2cTofNSigmaOuter, pi2cTofNSigmaOuter, float);
DECLARE_SOA_COLUMN(Pi2cHasRichSignal, pi2cHasRichSignal, bool);
DECLARE_SOA_COLUMN(Pi2cRichNSigma, pi2cRichNSigma, float);
DECLARE_SOA_COLUMN(Pi2cPdgCode, pi2cPdgCode, int);

DECLARE_SOA_COLUMN(PiccTofDeltaInner, piccTofDeltaInner, float);
DECLARE_SOA_COLUMN(PiccTofNSigmaInner, piccTofNSigmaInner, float);
DECLARE_SOA_COLUMN(PiccTofDeltaOuter, piccTofDeltaOuter, float);
DECLARE_SOA_COLUMN(PiccTofNSigmaOuter, piccTofNSigmaOuter, float);
DECLARE_SOA_COLUMN(PiccHasRichSignal, piccHasRichSignal, bool);
DECLARE_SOA_COLUMN(PiccRichNSigma, piccRichNSigma, float);
DECLARE_SOA_COLUMN(PiccPdgCode, piccPdgCode, int);

// Daughter info
DECLARE_SOA_COLUMN(PosPt, posPt, float);
DECLARE_SOA_COLUMN(PosEta, posEta, float);
DECLARE_SOA_COLUMN(NegPt, negPt, float);
DECLARE_SOA_COLUMN(NegEta, negEta, float);
DECLARE_SOA_COLUMN(BachPt, bachPt, float);
DECLARE_SOA_COLUMN(BachEta, bachEta, float);
DECLARE_SOA_COLUMN(BachPhi, bachPhi, float);
DECLARE_SOA_COLUMN(Pi1cPt, pi1cPt, float);
DECLARE_SOA_COLUMN(Pi1cEta, pi1cEta, float);
DECLARE_SOA_COLUMN(Pi2cPt, pi2cPt, float);
DECLARE_SOA_COLUMN(Pi2cEta, pi2cEta, float);
DECLARE_SOA_COLUMN(PiccPt, piccPt, float);
DECLARE_SOA_COLUMN(PiccEta, piccEta, float);

} // namespace otfmulticharm

DECLARE_SOA_TABLE(MCharmIndices, "AOD", "MCharmIndices",
                  o2::soa::Index<>,
                  otfmulticharm::CascadeId,
                  otfmulticharm::XiCPion1Id,
                  otfmulticharm::XiCPion2Id,
                  otfmulticharm::XiCCPionId);

DECLARE_SOA_TABLE(MCharmCores, "AOD", "MCharmCores",
                  otfmulticharm::XiccMass,
                  otfmulticharm::XiccPt,
                  otfmulticharm::XiccEta,
                  otfmulticharm::XiccDauDCA,

                  otfmulticharm::XicMass,
                  otfmulticharm::XicPt,
                  otfmulticharm::XicEta,
                  otfmulticharm::XicDauDCA,

                  otfmulticharm::XiDCAxy,
                  otfmulticharm::XiDCAz,
                  otfmulticharm::XicDCAxy,
                  otfmulticharm::XicDCAz,
                  otfmulticharm::XiccDCAxy,
                  otfmulticharm::XiccDCAz,

                  otfmulticharm::Pi1cDCAxy,
                  otfmulticharm::Pi1cDCAz,
                  otfmulticharm::Pi2cDCAxy,
                  otfmulticharm::Pi2cDCAz,
                  otfmulticharm::PiccDCAxy,
                  otfmulticharm::PiccDCAz,

                  otfmulticharm::XicDecayRadius2D,
                  otfmulticharm::XiccDecayRadius2D,
                  otfmulticharm::XicProperLength,
                  otfmulticharm::XicDistanceFromPV,
                  otfmulticharm::XiccProperLength,
                  otfmulticharm::Pi1cPt,
                  otfmulticharm::Pi2cPt,
                  otfmulticharm::PiccPt,
                  otfmulticharm::LUTConfigId);

DECLARE_SOA_TABLE(MCharmPID, "AOD", "MCharmPID",
                  otfmulticharm::Pi1cTofDeltaInner,
                  otfmulticharm::Pi1cTofNSigmaInner,
                  otfmulticharm::Pi1cTofDeltaOuter,
                  otfmulticharm::Pi1cTofNSigmaOuter,
                  otfmulticharm::Pi1cHasRichSignal,
                  otfmulticharm::Pi1cRichNSigma,
                  otfmulticharm::Pi1cPdgCode,

                  otfmulticharm::Pi2cTofDeltaInner,
                  otfmulticharm::Pi2cTofNSigmaInner,
                  otfmulticharm::Pi2cTofDeltaOuter,
                  otfmulticharm::Pi2cTofNSigmaOuter,
                  otfmulticharm::Pi2cHasRichSignal,
                  otfmulticharm::Pi2cRichNSigma,
                  otfmulticharm::Pi2cPdgCode,

                  otfmulticharm::PiccTofDeltaInner,
                  otfmulticharm::PiccTofNSigmaInner,
                  otfmulticharm::PiccTofDeltaOuter,
                  otfmulticharm::PiccTofNSigmaOuter,
                  otfmulticharm::PiccHasRichSignal,
                  otfmulticharm::PiccRichNSigma,
                  otfmulticharm::PiccPdgCode);

DECLARE_SOA_TABLE(MCharmExtra, "AOD", "MCharmExtra",
                  otfmulticharm::BachPt,
                  otfmulticharm::BachEta,
                  otfmulticharm::BachDCAxy,
                  otfmulticharm::BachDCAz,

                  otfmulticharm::PosPt,
                  otfmulticharm::PosEta,
                  otfmulticharm::PosDCAxy,
                  otfmulticharm::PosDCAz,

                  otfmulticharm::NegPt,
                  otfmulticharm::NegEta,
                  otfmulticharm::NegDCAxy,
                  otfmulticharm::NegDCAz,

                  otfmulticharm::Pi1cEta,
                  otfmulticharm::Pi2cEta,
                  otfmulticharm::PiccEta);

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFMULTICHARM_H_
