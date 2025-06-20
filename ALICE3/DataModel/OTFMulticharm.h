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

// kine vars
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);

// topo vars
DECLARE_SOA_COLUMN(XiDCAz, xiDCAz, float);
DECLARE_SOA_COLUMN(XiDCAxy, xiDCAxy, float);
DECLARE_SOA_COLUMN(XicDauDCA, xicDauDCA, float);
DECLARE_SOA_COLUMN(XicDCAxy, xicDCAxy, float);
DECLARE_SOA_COLUMN(XicDCAz, xicDCAz, float);
DECLARE_SOA_COLUMN(XiccDauDCA, xiccDauDCA, float);
DECLARE_SOA_COLUMN(XiccDCAxy, xiccDCAxy, float);
DECLARE_SOA_COLUMN(XiccDCAz, xiccDCAz, float);

DECLARE_SOA_COLUMN(PiFromXiDCAxy, piFromXiDCAxy, float);
DECLARE_SOA_COLUMN(PiFromLaDCAxy, piFromLaDCAxy, float);
DECLARE_SOA_COLUMN(PrFromLaDCAxy, prFromLaDCAxy, float);
DECLARE_SOA_COLUMN(PiFromXiDCAz, piFromXiDCAz, float);
DECLARE_SOA_COLUMN(PiFromLaDCAz, piFromLaDCAz, float);
DECLARE_SOA_COLUMN(PrFromLaDCAz, prFromLaDCAz, float);

DECLARE_SOA_COLUMN(Pi1cDCAxy, pi1cDCAxy, float);
DECLARE_SOA_COLUMN(Pi2cDCAxy, pi2cDCAxy, float);
DECLARE_SOA_COLUMN(PiccDCAxy, piccDCAxy, float);
DECLARE_SOA_COLUMN(Pi1cDCAz, pi1cDCAz, float);
DECLARE_SOA_COLUMN(Pi2cDCAz, pi2cDCAz, float);
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
DECLARE_SOA_COLUMN(Pi2cTofDeltaInner, pi2cTofDeltaInner, float);
DECLARE_SOA_COLUMN(Pi2cTofNSigmaInner, pi2cTofNSigmaInner, float);
DECLARE_SOA_COLUMN(Pi2cTofDeltaOuter, pi2cTofDeltaOuter, float);
DECLARE_SOA_COLUMN(Pi2cTofNSigmaOuter, pi2cTofNSigmaOuter, float);
DECLARE_SOA_COLUMN(PiccTofDeltaInner, piccTofDeltaInner, float);
DECLARE_SOA_COLUMN(PiccTofNSigmaInner, piccTofNSigmaInner, float);
DECLARE_SOA_COLUMN(PiccTofDeltaOuter, piccTofDeltaOuter, float);
DECLARE_SOA_COLUMN(PiccTofNSigmaOuter, piccTofNSigmaOuter, float);

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
                  otfmulticharm::XicDauDCA,
                  otfmulticharm::XiccDauDCA,
                  otfmulticharm::XicMass,
                  otfmulticharm::XiccMass,
                  otfmulticharm::Pt,
                  otfmulticharm::Eta,

                  otfmulticharm::XiDCAxy,
                  otfmulticharm::XiDCAz,
                  otfmulticharm::XicDCAxy,
                  otfmulticharm::XicDCAz,
                  otfmulticharm::XiccDCAxy,
                  otfmulticharm::XiccDCAz,

                  otfmulticharm::PiFromXiDCAxy,
                  otfmulticharm::PiFromXiDCAz,
                  otfmulticharm::PiFromLaDCAxy,
                  otfmulticharm::PiFromLaDCAz,
                  otfmulticharm::PrFromLaDCAxy,
                  otfmulticharm::PrFromLaDCAz,

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

                  otfmulticharm::Pi1cTofDeltaInner,
                  otfmulticharm::Pi1cTofNSigmaInner,
                  otfmulticharm::Pi1cTofDeltaOuter,
                  otfmulticharm::Pi1cTofNSigmaOuter,

                  otfmulticharm::Pi2cTofDeltaInner,
                  otfmulticharm::Pi2cTofNSigmaInner,
                  otfmulticharm::Pi2cTofDeltaOuter,
                  otfmulticharm::Pi2cTofNSigmaOuter,

                  otfmulticharm::PiccTofDeltaInner,
                  otfmulticharm::PiccTofNSigmaInner,
                  otfmulticharm::PiccTofDeltaOuter,
                  otfmulticharm::PiccTofNSigmaOuter,

                  otfmulticharm::BachPt,
                  otfmulticharm::BachEta,

                  otfmulticharm::PosPt,
                  otfmulticharm::PosEta,

                  otfmulticharm::NegPt,
                  otfmulticharm::NegEta,

                  otfmulticharm::Pi1cPt,
                  otfmulticharm::Pi1cEta,

                  otfmulticharm::Pi2cPt,
                  otfmulticharm::Pi2cEta,

                  otfmulticharm::PiccPt,
                  otfmulticharm::PiccEta);

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFMULTICHARM_H_
