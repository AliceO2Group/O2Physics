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
/// \file LFNonPromptCascadeTable.h
/// \brief Non prompt cascade tables
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_
#define PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_

namespace o2::aod
{
namespace NonPromptCascadeTableNS /// Per cosa sta NS? // come si aggiungono nuove variabili? (c'è una qualche lista?)
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);
DECLARE_SOA_COLUMN(Beta, beta, float);
DECLARE_SOA_COLUMN(Zvertex, zVertex, float);
DECLARE_SOA_COLUMN(DCAxy, dcaxy, float);
DECLARE_SOA_COLUMN(DCAz, dcaz, float);
DECLARE_SOA_COLUMN(TPCsignal, tpcSignal, float);
DECLARE_SOA_COLUMN(ITSchi2, itsChi2, float);
DECLARE_SOA_COLUMN(TPCchi2, tpcChi2, float);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(TPCfindableCls, tpcFindableCls, uint8_t);
DECLARE_SOA_COLUMN(TPCcrossedRows, tpcCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(ITSclsMap, itsClsMap, uint8_t);
DECLARE_SOA_COLUMN(TPCnCls, tpcNCls, uint8_t);
DECLARE_SOA_COLUMN(gPt, genPt, float);
DECLARE_SOA_COLUMN(gEta, genEta, float);
DECLARE_SOA_COLUMN(gPhi, genPhi, float);
DECLARE_SOA_COLUMN(PDGcode, pdgCode, int);
} // namespace NonPromptCascadeTableNS

DECLARE_SOA_TABLE(NonPromptCascadeTable, "AOD", "NONPROMPTCASCADETABLE",
                  NonPromptCascadeTableNS::Pt,
                  NonPromptCascadeTableNS::Eta,
                  NonPromptCascadeTableNS::Phi,
                  NonPromptCascadeTableNS::TPCInnerParam,
                  NonPromptCascadeTableNS::Beta,
                  NonPromptCascadeTableNS::Zvertex,
                  NonPromptCascadeTableNS::DCAxy,
                  NonPromptCascadeTableNS::DCAz,
                  NonPromptCascadeTableNS::TPCsignal,
                  NonPromptCascadeTableNS::ITSchi2,
                  NonPromptCascadeTableNS::TPCchi2,
                  NonPromptCascadeTableNS::Flags,
                  NonPromptCascadeTableNS::TPCfindableCls,
                  NonPromptCascadeTableNS::TPCcrossedRows,
                  NonPromptCascadeTableNS::ITSclsMap,
                  NonPromptCascadeTableNS::TPCnCls)

DECLARE_SOA_TABLE(NonPromptCascadeTableMC, "AOD", "NONPROMPTCASCADETABLEMC",
                  NonPromptCascadeTableNS::Pt,
                  NonPromptCascadeTableNS::Eta,
                  NonPromptCascadeTableNS::Phi,
                  NonPromptCascadeTableNS::TPCInnerParam,
                  NonPromptCascadeTableNS::Beta,
                  NonPromptCascadeTableNS::Zvertex,
                  NonPromptCascadeTableNS::DCAxy,
                  NonPromptCascadeTableNS::DCAz,
                  NonPromptCascadeTableNS::TPCsignal,
                  NonPromptCascadeTableNS::ITSchi2,
                  NonPromptCascadeTableNS::TPCchi2,
                  NonPromptCascadeTableNS::Flags,
                  NonPromptCascadeTableNS::TPCfindableCls,
                  NonPromptCascadeTableNS::TPCcrossedRows,
                  NonPromptCascadeTableNS::ITSclsMap,
                  NonPromptCascadeTableNS::TPCnCls,
                  NonPromptCascadeTableNS::gPt,
                  NonPromptCascadeTableNS::gEta,
                  NonPromptCascadeTableNS::gPhi,
                  NonPromptCascadeTableNS::PDGcode)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_
