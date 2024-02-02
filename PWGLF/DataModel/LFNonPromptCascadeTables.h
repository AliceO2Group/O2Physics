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
namespace NPCascadeTable
{
DECLARE_SOA_COLUMN(CascPt, cascPt, float);
DECLARE_SOA_COLUMN(CascEta, cascEta, float);
DECLARE_SOA_COLUMN(CascPhi, cascPhi, float);

DECLARE_SOA_COLUMN(CascDCAxy, cascDCAxy, float);
DECLARE_SOA_COLUMN(CascDCAz, cascDCAz, float);
DECLARE_SOA_COLUMN(ProtonDCAxy, protonDCAxy, float);
DECLARE_SOA_COLUMN(ProtonDCAz, protonDCAz, float);
DECLARE_SOA_COLUMN(PionDCAxy, pionDCAxy, float);
DECLARE_SOA_COLUMN(PionDCAz, pionDCAz, float);
DECLARE_SOA_COLUMN(BachDCAxy, bachDCAxy, float);
DECLARE_SOA_COLUMN(BachDCAz, bachDCAz, float);

DECLARE_SOA_COLUMN(CascCosPA, casccosPA, float);
DECLARE_SOA_COLUMN(V0CosPA, v0cosPA, float);

DECLARE_SOA_COLUMN(MassXi, massXi, double);
DECLARE_SOA_COLUMN(MassOmega, massOmega, double);
DECLARE_SOA_COLUMN(MassV0, massV0, double);

DECLARE_SOA_COLUMN(CascRadius, cascRadius, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);

DECLARE_SOA_COLUMN(CascLenght, cascLenght, float);
DECLARE_SOA_COLUMN(V0Lenght, v0Lenght, float);

DECLARE_SOA_COLUMN(CascNClusITS, cascNClusITS, int);
DECLARE_SOA_COLUMN(ProtonNClusITS, protonNClusITS, int);
DECLARE_SOA_COLUMN(PionNClusITS, pionNClusITS, int);
DECLARE_SOA_COLUMN(BachKaonNClusITS, bachKaonNClusITS, int);
DECLARE_SOA_COLUMN(BachPionNClusITS, bachPionNClusITS, int);

DECLARE_SOA_COLUMN(ProtonNClusTPC, protonNClusTPC, int);
DECLARE_SOA_COLUMN(PionNClusTPC, pionNClusTPC, int);
DECLARE_SOA_COLUMN(BachKaonNClusTPC, bachKaonNClusTPC, int);
DECLARE_SOA_COLUMN(BachPionNClusTPC, bachPionNClusTPC, int);

DECLARE_SOA_COLUMN(ProtonTPCNSigma, protonTPCNSigma, float);
DECLARE_SOA_COLUMN(PionTPCNSigma, pionTPCNSigma, float);
DECLARE_SOA_COLUMN(BachKaonTPCNSigma, bachKaonTPCNSigma, float);
DECLARE_SOA_COLUMN(BachPionTPCNSigma, bachPionTPCNSigma, float);

DECLARE_SOA_COLUMN(ProtonHasTOF, protonHasTOF, bool);
DECLARE_SOA_COLUMN(PionHasTOF, pionHasTOF, bool);
DECLARE_SOA_COLUMN(BachKaonHasTOF, bachKaonHasTOF, bool);
DECLARE_SOA_COLUMN(BachPionHasTOF, bachPionHasTOF, bool);

DECLARE_SOA_COLUMN(ProtonTOFNSigma, protonTOFNSigma, float);
DECLARE_SOA_COLUMN(PionTOFNSigma, pionTOFNSigma, float);
DECLARE_SOA_COLUMN(BachKaonTOFNSigma, bachKaonTOFNSigma, float);
DECLARE_SOA_COLUMN(BachPionTOFNSigma, bachPionTOFNSigma, float);

DECLARE_SOA_COLUMN(gPt, genPt, float);
DECLARE_SOA_COLUMN(gEta, genEta, float);
DECLARE_SOA_COLUMN(gPhi, genPhi, float);
DECLARE_SOA_COLUMN(PDGcode, pdgCode, int);

} // namespace NPCascadeTable
DECLARE_SOA_TABLE(NPCascTable, "AOD", "NPCASCTABLE",
                  NPCascadeTable::CascPt,
                  NPCascadeTable::CascEta,
                  NPCascadeTable::CascPhi,
                  NPCascadeTable::CascDCAxy,
                  NPCascadeTable::CascDCAz,
                  NPCascadeTable::ProtonDCAxy,
                  NPCascadeTable::ProtonDCAz,
                  NPCascadeTable::PionDCAxy,
                  NPCascadeTable::PionDCAz,
                  NPCascadeTable::BachDCAxy,
                  NPCascadeTable::BachDCAz,
                  NPCascadeTable::CascCosPA,
                  NPCascadeTable::V0CosPA,
                  NPCascadeTable::MassXi,
                  NPCascadeTable::MassOmega,
                  NPCascadeTable::MassV0,
                  NPCascadeTable::CascRadius,
                  NPCascadeTable::V0Radius,
                  NPCascadeTable::CascLenght,
                  NPCascadeTable::V0Lenght,
                  NPCascadeTable::CascNClusITS,
                  NPCascadeTable::ProtonNClusITS,
                  NPCascadeTable::PionNClusITS,
                  NPCascadeTable::BachKaonNClusITS,
                  NPCascadeTable::BachPionNClusITS,
                  NPCascadeTable::ProtonNClusTPC,
                  NPCascadeTable::PionNClusTPC,
                  NPCascadeTable::BachKaonNClusTPC,
                  NPCascadeTable::BachPionNClusTPC,
                  NPCascadeTable::ProtonTPCNSigma,
                  NPCascadeTable::PionTPCNSigma,
                  NPCascadeTable::BachKaonTPCNSigma,
                  NPCascadeTable::BachPionTPCNSigma,
                  NPCascadeTable::ProtonHasTOF,
                  NPCascadeTable::PionHasTOF,
                  NPCascadeTable::BachKaonHasTOF,
                  NPCascadeTable::BachPionHasTOF,
                  NPCascadeTable::ProtonTOFNSigma,
                  NPCascadeTable::PionTOFNSigma,
                  NPCascadeTable::BachKaonTOFNSigma,
                  NPCascadeTable::BachPionTOFNSigma)

DECLARE_SOA_TABLE(NPCascTableMC, "AOD", "NPCASCTABLEMC",
                  NPCascadeTable::CascPt,
                  NPCascadeTable::CascEta,
                  NPCascadeTable::CascPhi,
                  NPCascadeTable::CascDCAxy,
                  NPCascadeTable::CascDCAz,
                  NPCascadeTable::ProtonDCAxy,
                  NPCascadeTable::ProtonDCAz,
                  NPCascadeTable::PionDCAxy,
                  NPCascadeTable::PionDCAz,
                  NPCascadeTable::BachDCAxy,
                  NPCascadeTable::BachDCAz,
                  NPCascadeTable::CascCosPA,
                  NPCascadeTable::V0CosPA,
                  NPCascadeTable::MassXi,
                  NPCascadeTable::MassOmega,
                  NPCascadeTable::MassV0,
                  NPCascadeTable::CascRadius,
                  NPCascadeTable::V0Radius,
                  NPCascadeTable::CascLenght,
                  NPCascadeTable::V0Lenght,
                  NPCascadeTable::CascNClusITS,
                  NPCascadeTable::ProtonNClusITS,
                  NPCascadeTable::PionNClusITS,
                  NPCascadeTable::BachKaonNClusITS,
                  NPCascadeTable::BachPionNClusITS,
                  NPCascadeTable::ProtonNClusTPC,
                  NPCascadeTable::PionNClusTPC,
                  NPCascadeTable::BachKaonNClusTPC,
                  NPCascadeTable::BachPionNClusTPC,
                  NPCascadeTable::ProtonTPCNSigma,
                  NPCascadeTable::PionTPCNSigma,
                  NPCascadeTable::BachKaonTPCNSigma,
                  NPCascadeTable::BachPionTPCNSigma,
                  NPCascadeTable::ProtonHasTOF,
                  NPCascadeTable::PionHasTOF,
                  NPCascadeTable::BachKaonHasTOF,
                  NPCascadeTable::BachPionHasTOF,
                  NPCascadeTable::ProtonTOFNSigma,
                  NPCascadeTable::PionTOFNSigma,
                  NPCascadeTable::BachKaonTOFNSigma,
                  NPCascadeTable::BachPionTOFNSigma,
                  NPCascadeTable::gPt,
                  NPCascadeTable::gEta,
                  NPCascadeTable::gPhi,
                  NPCascadeTable::PDGcode)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_
