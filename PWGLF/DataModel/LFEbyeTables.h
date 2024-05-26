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

namespace LFEbyeTable
{
DECLARE_SOA_COLUMN(IdColl, idColl, int64_t);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Zvtx, zvtx, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);
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
DECLARE_SOA_COLUMN(GenPt, gentPt, float);
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(IsReco, isReco, bool);
} // namespace LFEbyeTable

DECLARE_SOA_TABLE(NucleiEbyeTable, "AOD", "NUCLEBYETABLE",
                  LFEbyeTable::IdColl,
                  LFEbyeTable::Centrality,
                  LFEbyeTable::Zvtx,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::IsMatter,
                  LFEbyeTable::DcaPV,
                  LFEbyeTable::TpcNcls,
                  LFEbyeTable::TpcNsigma,
                  LFEbyeTable::TofMass);

DECLARE_SOA_TABLE(McNucleiEbyeTable, "AOD", "MCNUCLEBYETABLE",
                  LFEbyeTable::IdColl,
                  LFEbyeTable::Centrality,
                  LFEbyeTable::Zvtx,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::IsMatter,
                  LFEbyeTable::DcaPV,
                  LFEbyeTable::TpcNcls,
                  LFEbyeTable::TpcNsigma,
                  LFEbyeTable::TofMass,
                  LFEbyeTable::GenPt,
                  LFEbyeTable::GenEta,
                  LFEbyeTable::PdgCode,
                  LFEbyeTable::IsReco);

DECLARE_SOA_TABLE(LambdaEbyeTable, "AOD", "LAMBEBYETABLE",
                  LFEbyeTable::IdColl,
                  LFEbyeTable::Centrality,
                  LFEbyeTable::Zvtx,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::IsMatter,
                  LFEbyeTable::DcaV0PV,
                  LFEbyeTable::DcaNegPV,
                  LFEbyeTable::DcaPosPV,
                  LFEbyeTable::DcaV0Tracks,
                  LFEbyeTable::CosPA,
                  LFEbyeTable::TpcNsigmaNeg,
                  LFEbyeTable::TpcNsigmaPos,
                  LFEbyeTable::IdNeg,
                  LFEbyeTable::IdPos);

DECLARE_SOA_TABLE(McLambdaEbyeTable, "AOD", "MCLAMBEBYETABLE",
                  LFEbyeTable::IdColl,
                  LFEbyeTable::Centrality,
                  LFEbyeTable::Zvtx,
                  LFEbyeTable::Pt,
                  LFEbyeTable::Eta,
                  LFEbyeTable::Mass,
                  LFEbyeTable::IsMatter,
                  LFEbyeTable::DcaV0PV,
                  LFEbyeTable::DcaNegPV,
                  LFEbyeTable::DcaPosPV,
                  LFEbyeTable::DcaV0Tracks,
                  LFEbyeTable::CosPA,
                  LFEbyeTable::TpcNsigmaNeg,
                  LFEbyeTable::TpcNsigmaPos,
                  LFEbyeTable::IdNeg,
                  LFEbyeTable::IdPos,
                  LFEbyeTable::GenPt,
                  LFEbyeTable::GenEta,
                  LFEbyeTable::PdgCode,
                  LFEbyeTable::IsReco);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFEBYENUCLEISTRANGE_H_
