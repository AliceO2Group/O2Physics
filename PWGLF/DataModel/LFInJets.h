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
/// \brief Derived Data table for LF in jets analysis
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \author Sara Pucillo (sara.pucillo@cern.ch)

#ifndef PWGLF_DATAMODEL_LFINJETS_H_
#define PWGLF_DATAMODEL_LFINJETS_H_

#include <Framework/ASoA.h>

namespace o2::aod
{

namespace lfinjets
{

DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Sign, sign, int);
DECLARE_SOA_COLUMN(MassLambda, masslambda, float);
DECLARE_SOA_COLUMN(MassAntiLambda, massantilambda, float);
DECLARE_SOA_COLUMN(MassK0Short, massk0short, float);
DECLARE_SOA_COLUMN(V0Radius, v0radius, float);
DECLARE_SOA_COLUMN(V0CosPA, v0cospa, float);
DECLARE_SOA_COLUMN(V0DCAPosToPV, v0dcapostopv, float);
DECLARE_SOA_COLUMN(V0DCANegToPV, v0dcanegtopv, float);
DECLARE_SOA_COLUMN(V0DCAV0Daughters, v0dcav0daughters, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPr, ntpcsigmanegpr, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPr, ntpcsigmapospr, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPi, ntpcsigmanegpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPi, ntpcsigmapospi, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPr, ntofsigmanegpr, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPr, ntofsigmapospr, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPi, ntofsigmanegpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPi, ntofsigmapospi, float);
DECLARE_SOA_COLUMN(MultFT0M, multft0m, float);
DECLARE_SOA_COLUMN(V0PosTPCCrossedRows, v0postpcCrossedRows, float);
DECLARE_SOA_COLUMN(V0NegTPCCrossedRows, v0negtpcCrossedRows, float);
DECLARE_SOA_COLUMN(V0NegTPCChi2, v0negTPCChi2, float);
DECLARE_SOA_COLUMN(V0NegITSlayers, v0negITSlayers, int);
DECLARE_SOA_COLUMN(V0PosTPCChi2, v0posTPCChi2, float);
DECLARE_SOA_COLUMN(V0PosITSlayers, v0posITSlayers, int);
DECLARE_SOA_COLUMN(MassXi, massxi, float);
DECLARE_SOA_COLUMN(MassOmega, massomega, float);
DECLARE_SOA_COLUMN(CascRadius, cascradius, float);
DECLARE_SOA_COLUMN(CascCosPA, casccospa, float);
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);
DECLARE_SOA_COLUMN(DCACascDaughters, dcacascdaughters, float);
DECLARE_SOA_COLUMN(DCAV0ToPV, dcav0topv, float);
DECLARE_SOA_COLUMN(NTPCSigmaBachPi, ntpcsigmabachpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaBachKa, ntpcsigmabachka, float);
DECLARE_SOA_COLUMN(NTOFSigmaBachPi, ntofsigmabachpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaBachKa, ntofsigmabachka, float);
DECLARE_SOA_COLUMN(BachTPCCrossedRows, bachtpcCrossedRows, float);
DECLARE_SOA_COLUMN(BachTPCChi2, bachTPCChi2, float);
DECLARE_SOA_COLUMN(BachITSlayers, bachITSlayers, int);
DECLARE_SOA_COLUMN(IsUE, isUE, bool);
DECLARE_SOA_COLUMN(IsJC, isJC, bool);

} // namespace lfinjets

DECLARE_SOA_TABLE(V0InJets, "AOD", "V0INJETS",
                  lfinjets::Pt, lfinjets::MassLambda, lfinjets::MassAntiLambda, lfinjets::MassK0Short,
                  lfinjets::V0Radius, lfinjets::V0CosPA, lfinjets::V0DCAPosToPV,
                  lfinjets::V0DCANegToPV, lfinjets::V0DCAV0Daughters,
                  lfinjets::NTPCSigmaNegPr, lfinjets::NTPCSigmaPosPr, lfinjets::NTPCSigmaNegPi, lfinjets::NTPCSigmaPosPi,
                  lfinjets::NTOFSigmaNegPr, lfinjets::NTOFSigmaPosPr, lfinjets::NTOFSigmaNegPi, lfinjets::NTOFSigmaPosPi,
                  lfinjets::MultFT0M, lfinjets::V0PosTPCCrossedRows, lfinjets::V0NegTPCCrossedRows,
                  lfinjets::V0NegTPCChi2, lfinjets::V0NegITSlayers, lfinjets::V0PosTPCChi2, lfinjets::V0PosITSlayers,
                  lfinjets::IsUE, lfinjets::IsJC);

DECLARE_SOA_TABLE(CascInJets, "AOD", "CASCINJETS",
                  lfinjets::Pt, lfinjets::Sign, lfinjets::MassXi, lfinjets::MassOmega, lfinjets::MassLambda,
                  lfinjets::CascRadius, lfinjets::CascCosPA, lfinjets::V0Radius, lfinjets::V0CosPA,
                  lfinjets::V0DCAPosToPV, lfinjets::V0DCANegToPV, lfinjets::DCABachToPV,
                  lfinjets::DCACascDaughters, lfinjets::V0DCAV0Daughters, lfinjets::DCAV0ToPV,
                  lfinjets::NTPCSigmaNegPr, lfinjets::NTPCSigmaPosPr, lfinjets::NTPCSigmaNegPi, lfinjets::NTPCSigmaPosPi,
                  lfinjets::NTPCSigmaBachPi, lfinjets::NTPCSigmaBachKa,
                  lfinjets::NTOFSigmaNegPr, lfinjets::NTOFSigmaPosPr, lfinjets::NTOFSigmaNegPi, lfinjets::NTOFSigmaPosPi,
                  lfinjets::NTOFSigmaBachPi, lfinjets::NTOFSigmaBachKa, lfinjets::MultFT0M,
                  lfinjets::V0PosTPCCrossedRows, lfinjets::V0NegTPCCrossedRows, lfinjets::BachTPCCrossedRows,
                  lfinjets::V0NegTPCChi2, lfinjets::V0NegITSlayers, lfinjets::V0PosTPCChi2, lfinjets::V0PosITSlayers,
                  lfinjets::BachTPCChi2, lfinjets::BachITSlayers,
                  lfinjets::IsUE, lfinjets::IsJC);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFINJETS_H_
