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
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFEBYETABLES_H_
#define PWGLF_DATAMODEL_LFEBYETABLES_H_

namespace o2::aod
{

namespace LFEbyeCollTable
{
DECLARE_SOA_COLUMN(Centrality, centrality, uint8_t);
DECLARE_SOA_COLUMN(Zvtx, zvtx, float);
DECLARE_SOA_COLUMN(ZvtxMask, zvtxMask, int8_t);
DECLARE_SOA_COLUMN(TriggerMask, triggerMask, uint8_t);
DECLARE_SOA_COLUMN(CBMultiplicity, cbMultiplicity, uint8_t);
DECLARE_SOA_COLUMN(Ntracks, ntracks, uint8_t);
} // namespace LFEbyeCollTable

DECLARE_SOA_TABLE(CollEbyeTables, "AOD", "COLLEBYETABLE",
                  o2::soa::Index<>,
                  LFEbyeCollTable::Centrality,
                  LFEbyeCollTable::Zvtx);
using CollEbyeTable = CollEbyeTables::iterator;

DECLARE_SOA_TABLE(MiniCollTables, "AOD", "MINICOLLTABLE",
                  o2::soa::Index<>,
                  LFEbyeCollTable::ZvtxMask,
                  LFEbyeCollTable::TriggerMask,
                  LFEbyeCollTable::CBMultiplicity,
                  LFEbyeCollTable::Centrality,
                  LFEbyeCollTable::Ntracks);
using MiniCollTable = MiniCollTables::iterator;

namespace LFEbyeTable
{
DECLARE_SOA_INDEX_COLUMN(CollEbyeTable, collEbyeTable);
DECLARE_SOA_INDEX_COLUMN(MiniCollTable, miniCollTable);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(DcaPV, dcaPV, float);
DECLARE_SOA_COLUMN(TpcNcls, tpcNcls, uint8_t);
DECLARE_SOA_COLUMN(TpcNsigma, tpcNsigma, float);
DECLARE_SOA_COLUMN(TofMass, tofMass, float);
DECLARE_SOA_COLUMN(DcaV0PV, dcaV0Pv, float);
DECLARE_SOA_COLUMN(DcaNegPV, dcaNegPv, float);
DECLARE_SOA_COLUMN(DcaPosPV, dcaPosPv, float);
DECLARE_SOA_COLUMN(DcaV0Tracks, dcaV0tracks, float);
DECLARE_SOA_COLUMN(CosPA, cosPa, double);
DECLARE_SOA_COLUMN(TpcNsigmaNeg, tpcNsigmaNeg, float);
DECLARE_SOA_COLUMN(TpcNsigmaPos, tpcNsigmaPos, float);
DECLARE_SOA_COLUMN(IdNeg, idNeg, int64_t);
DECLARE_SOA_COLUMN(IdPos, idPos, int64_t);
DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(IsReco, isReco, bool);
DECLARE_SOA_COLUMN(EtaMask, etaMask, int8_t);
DECLARE_SOA_COLUMN(SelMask, selMask, int);
DECLARE_SOA_COLUMN(OuterPID, outerPID, float);
DECLARE_SOA_COLUMN(GenEtaMask, genEtaMask, int8_t);
} // namespace LFEbyeTable

DECLARE_SOA_TABLE(NucleiEbyeTables, "AOD", "NUCLEBYETABLE",
                  o2::soa::Index<>,
                  LFEbyeTable::CollEbyeTableId,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::DcaPV,
                  LFEbyeTable::TpcNcls,
                  LFEbyeTable::TpcNsigma,
                  LFEbyeTable::TofMass);
using NucleiEbyeTable = NucleiEbyeTables::iterator;

DECLARE_SOA_TABLE(McNucleiEbyeTables, "AOD", "MCNUCLEBYETABLE",
                  o2::soa::Index<>,
                  LFEbyeTable::CollEbyeTableId,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::DcaPV,
                  LFEbyeTable::TpcNcls,
                  LFEbyeTable::TpcNsigma,
                  LFEbyeTable::TofMass,
                  LFEbyeTable::GenPt,
                  LFEbyeTable::GenEta,
                  LFEbyeTable::PdgCode,
                  LFEbyeTable::IsReco);
using McNucleiEbyeTable = McNucleiEbyeTables::iterator;

DECLARE_SOA_TABLE(LambdaEbyeTables, "AOD", "LAMBEBYETABLE",
                  o2::soa::Index<>,
                  LFEbyeTable::CollEbyeTableId,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::DcaV0PV,
                  // LFEbyeTable::DcaNegPV,
                  // LFEbyeTable::DcaPosPV,
                  LFEbyeTable::DcaV0Tracks,
                  LFEbyeTable::CosPA,
                  // LFEbyeTable::TpcNsigmaNeg,
                  // LFEbyeTable::TpcNsigmaPos,
                  LFEbyeTable::IdNeg,
                  LFEbyeTable::IdPos);
using LambdaEbyeTable = LambdaEbyeTables::iterator;

DECLARE_SOA_TABLE(McLambdaEbyeTables, "AOD", "MCLAMBEBYETABLE",
                  o2::soa::Index<>,
                  LFEbyeTable::CollEbyeTableId,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::DcaV0PV,
                  // LFEbyeTable::DcaNegPV,
                  // LFEbyeTable::DcaPosPV,
                  LFEbyeTable::DcaV0Tracks,
                  LFEbyeTable::CosPA,
                  // LFEbyeTable::TpcNsigmaNeg,
                  // LFEbyeTable::TpcNsigmaPos,
                  LFEbyeTable::IdNeg,
                  LFEbyeTable::IdPos,
                  LFEbyeTable::GenPt,
                  LFEbyeTable::GenEta,
                  LFEbyeTable::PdgCode,
                  LFEbyeTable::IsReco);
using McLambdaEbyeTable = McLambdaEbyeTables::iterator;

DECLARE_SOA_TABLE(MiniTrkTables, "AOD", "MINITRKTABLE",
                  o2::soa::Index<>,
                  LFEbyeTable::MiniCollTableId,
                  LFEbyeTable::Pt,
                  LFEbyeTable::EtaMask,
                  LFEbyeTable::SelMask,
                  LFEbyeTable::OuterPID);
using MiniTrkTable = MiniTrkTables::iterator;

DECLARE_SOA_TABLE(McMiniTrkTables, "AOD", "MCMINITRKTABLE",
                  o2::soa::Index<>,
                  LFEbyeTable::MiniCollTableId,
                  LFEbyeTable::Pt,
                  LFEbyeTable::EtaMask,
                  LFEbyeTable::SelMask,
                  LFEbyeTable::OuterPID,
                  LFEbyeTable::GenPt,
                  LFEbyeTable::GenEtaMask,
                  LFEbyeTable::IsReco);
using McMiniTrkTable = McMiniTrkTables::iterator;
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFEBYETABLES_H_
