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
DECLARE_SOA_COLUMN(MatchingChi2, matchingChi2, float);
DECLARE_SOA_COLUMN(DeltaPtITSCascade, deltaPtITSCascade, float);
DECLARE_SOA_COLUMN(ITSClusSize, itsClusSize, float);
DECLARE_SOA_COLUMN(HasReassociatedCluster, hasReassociatedCluster, bool);
DECLARE_SOA_COLUMN(IsGoodMatch, isGoodMatch, bool);
DECLARE_SOA_COLUMN(IsGoodCascade, isGoodCascade, bool);
DECLARE_SOA_COLUMN(PdgCodeMom, pdgCodeMom, int);
DECLARE_SOA_COLUMN(PdgCodeITStrack, pdgCodeITStrack, int);
DECLARE_SOA_COLUMN(IsFromBeauty, isFromBeauty, bool);
DECLARE_SOA_COLUMN(IsFromCharm, isFromCharm, bool);

DECLARE_SOA_COLUMN(PvX, pvX, float);
DECLARE_SOA_COLUMN(PvY, pvY, float);
DECLARE_SOA_COLUMN(PvZ, pvZ, float);

DECLARE_SOA_COLUMN(CascPt, cascPt, float);
DECLARE_SOA_COLUMN(CascEta, cascEta, float);
DECLARE_SOA_COLUMN(CascPhi, cascPhi, float);

DECLARE_SOA_COLUMN(ProtonPt, protonPt, float);
DECLARE_SOA_COLUMN(ProtonEta, protonEta, float);
DECLARE_SOA_COLUMN(PionPt, pionPt, float);
DECLARE_SOA_COLUMN(PionEta, pionEta, float);
DECLARE_SOA_COLUMN(BachPt, bachPt, float);
DECLARE_SOA_COLUMN(BachEta, bachEta, float);

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
DECLARE_SOA_COLUMN(DCAxMC, dcaXmc, float);
DECLARE_SOA_COLUMN(DCAyMC, dcaYmc, float);
DECLARE_SOA_COLUMN(DCAzMC, dcaZmc, float);
DECLARE_SOA_COLUMN(MCcollisionMatch, mcCollisionMatch, bool);

} // namespace NPCascadeTable
DECLARE_SOA_TABLE(NPCascTable, "AOD", "NPCASCTABLE",
                  NPCascadeTable::MatchingChi2,
                  NPCascadeTable::DeltaPtITSCascade,
                  NPCascadeTable::ITSClusSize,
                  NPCascadeTable::HasReassociatedCluster,
                  aod::collision::NumContrib,
                  aod::collision::CollisionTimeRes,
                  NPCascadeTable::PvX,
                  NPCascadeTable::PvY,
                  NPCascadeTable::PvZ,
                  NPCascadeTable::CascPt,
                  NPCascadeTable::CascEta,
                  NPCascadeTable::CascPhi,
                  NPCascadeTable::ProtonPt,
                  NPCascadeTable::ProtonEta,
                  NPCascadeTable::PionPt,
                  NPCascadeTable::PionEta,
                  NPCascadeTable::BachPt,
                  NPCascadeTable::BachEta,
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
                  NPCascadeTable::MatchingChi2,
                  NPCascadeTable::DeltaPtITSCascade,
                  NPCascadeTable::ITSClusSize,
                  NPCascadeTable::HasReassociatedCluster,
                  NPCascadeTable::IsGoodMatch,
                  NPCascadeTable::IsGoodCascade,
                  NPCascadeTable::PdgCodeMom,
                  NPCascadeTable::PdgCodeITStrack,
                  NPCascadeTable::IsFromBeauty,
                  NPCascadeTable::IsFromCharm,
                  aod::collision::NumContrib,
                  aod::collision::CollisionTimeRes,
                  NPCascadeTable::PvX,
                  NPCascadeTable::PvY,
                  NPCascadeTable::PvZ,
                  NPCascadeTable::CascPt,
                  NPCascadeTable::CascEta,
                  NPCascadeTable::CascPhi,
                  NPCascadeTable::ProtonPt,
                  NPCascadeTable::ProtonEta,
                  NPCascadeTable::PionPt,
                  NPCascadeTable::PionEta,
                  NPCascadeTable::BachPt,
                  NPCascadeTable::BachEta,
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
                  NPCascadeTable::PDGcode,
                  NPCascadeTable::DCAxMC,
                  NPCascadeTable::DCAyMC,
                  NPCascadeTable::DCAzMC,
                  NPCascadeTable::MCcollisionMatch)

DECLARE_SOA_TABLE(NPCascTableGen, "AOD", "NPCASCTABLEGen",
                  NPCascadeTable::gPt,
                  NPCascadeTable::gEta,
                  NPCascadeTable::gPhi,
                  NPCascadeTable::PDGcode,
                  NPCascadeTable::PdgCodeMom,
                  NPCascadeTable::DCAxMC,
                  NPCascadeTable::DCAyMC,
                  NPCascadeTable::DCAzMC,
                  NPCascadeTable::IsFromBeauty,
                  NPCascadeTable::IsFromCharm)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_
