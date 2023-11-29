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
DECLARE_SOA_COLUMN(PionDCAxy,pionDCAxy, float);
DECLARE_SOA_COLUMN(PionDCAz, pionDCAz, float);
DECLARE_SOA_COLUMN(BachDCAxy, bachDCAxy, float);
DECLARE_SOA_COLUMN(BachDCAz, bachDCAz, float);

DECLARE_SOA_COLUMN(CascCosPA, casccosPA, float);
DECLARE_SOA_COLUMN(V0CosPA, v0cosPA, float);
DECLARE_SOA_COLUMN(BachCOSPA, bachcosPA, float);

DECLARE_SOA_COLUMN(MassXi, massXi, double);
DECLARE_SOA_COLUMN(MassOmega, massOmega, double);
DECLARE_SOA_COLUMN(MassV0, massV0, double);

DECLARE_SOA_COLUMN(CascRadius, cascRadius, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);

DECLARE_SOA_COLUMN(CascLenght, cascLenght, float);
DECLARE_SOA_COLUMN(V0Lenght, v0Lenght, float);

DECLARE_SOA_COLUMN(CascNClusITS, cascNClusITS, int);

DECLARE_SOA_COLUMN(CascNClusTPC, cascNClusTPC, int); 
DECLARE_SOA_COLUMN(ProtonNClusTPC, protonNClusTPC, int); 
DECLARE_SOA_COLUMN(PionNClusTPC, pionNClusTPC, int); 
DECLARE_SOA_COLUMN(BachKaonNClusTPC, bachKaonNClusTPC, int);
DECLARE_SOA_COLUMN(BachPionNClusTPC, bachPionNClusTPC, int);
DECLARE_SOA_COLUMN(ProtonTPCSignal, protonTPCSignal, float); 
DECLARE_SOA_COLUMN(PionTPCSignal, pionTPCSignal, float); 
DECLARE_SOA_COLUMN(BachKaonTPCSignal, bachKaonTPCSignal, float); 
DECLARE_SOA_COLUMN(BachPionTPCSignal, bachPionTPCSignal, float); 
DECLARE_SOA_COLUMN(ProtonTPCNSigma, protonTPCNSigma, float); 
DECLARE_SOA_COLUMN(PionTPCNSigma, pionTPCNSigma, float); 
DECLARE_SOA_COLUMN(BachKaonTPCNSigma, bachKaonTPCNSigma, float); 
DECLARE_SOA_COLUMN(BachPionTPCNSigma, bachPionTPCNSigma, float); 

DECLARE_SOA_COLUMN(ProtonHasTOF, protonHasTOF, bool); 
DECLARE_SOA_COLUMN(PionHasTOF, pionHasTOF, bool); 
DECLARE_SOA_COLUMN(BachKaonHasTOF, bachKaonHasTOF, bool); 
DECLARE_SOA_COLUMN(BachPionHasTOF, bachPionHasTOF, bool); 
DECLARE_SOA_COLUMN(ProtonTOFSignal, protonTOFSignal, float);
DECLARE_SOA_COLUMN(PionTOFSignal, pionTOFSignal, float);
DECLARE_SOA_COLUMN(BachKaonTOFSignal, bachKaonTOFSignal, float);
DECLARE_SOA_COLUMN(BachPionTOFSignal, bachPionTOFSignal, float);
DECLARE_SOA_COLUMN(ProtonTOFNSigma, protonTOFNSigma, float);
DECLARE_SOA_COLUMN(PionTOFNSigma, pionTOFNSigma, float);
DECLARE_SOA_COLUMN(BachKaonTOFNSigma, bachKaonTOFNSigma, float);
DECLARE_SOA_COLUMN(BachPionTOFNSigma, bachPionTOFNSigma, float);

DECLARE_SOA_COLUMN(gPt, genPt, float);
DECLARE_SOA_COLUMN(gEta, genEta, float);
DECLARE_SOA_COLUMN(gPhi, genPhi, float);
DECLARE_SOA_COLUMN(PDGcode, pdgCode, int);

} // namespace NPCascadeTable
DECLARE_SOA_TABLE(NonPromptCascadeTable, "AOD", "NPCASCTABLE",
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
                  NPCascadeTable::BachCOSPA,
                  NPCascadeTable::MassXi,
                  NPCascadeTable::MassOmega,
                  NPCascadeTable::MassV0,
                  NPCascadeTable::CascRadius,
                  NPCascadeTable::V0Radius,
                  NPCascadeTable::CascLenght,
                  NPCascadeTable::V0Lenght,
                  NPCascadeTable::CascNClusITS,
                  NPCascadeTable::CascNClusTPC,
                  NPCascadeTable::ProtonNClusTPC,
                  NPCascadeTable::PionNClusTPC,
                  NPCascadeTable::BachKaonNClusTPC,
                  NPCascadeTable::BachPionNClusTPC,
                  NPCascadeTable::ProtonTPCSignal,
                  NPCascadeTable::PionTPCSignal,
                  NPCascadeTable::BachKaonTPCSignal,
                  NPCascadeTable::BachPionTPCSignal,
                  NPCascadeTable::ProtonTPCNSigma,
                  NPCascadeTable::PionTPCNSigma,
                  NPCascadeTable::BachKaonTPCNSigma,
                  NPCascadeTable::BachPionTPCNSigma,
                  NPCascadeTable::ProtonHasTOF,
                  NPCascadeTable::PionHasTOF,
                  NPCascadeTable::BachKaonHasTOF,
                  NPCascadeTable::BachPionHasTOF,
                  NPCascadeTable::ProtonTOFSignal,
                  NPCascadeTable::PionTOFSignal,
                  NPCascadeTable::BachKaonTOFSignal,
                  NPCascadeTable::BachPionTOFSignal,
                  NPCascadeTable::ProtonTOFNSigma,
                  NPCascadeTable::PionTOFNSigma,
                  NPCascadeTable::BachKaonTOFNSigma,
                  NPCascadeTable::BachPionTOFNSigma)

DECLARE_SOA_TABLE(NonPromptCascadeTableMC, "AOD", "NPCASCTABLEMC",
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
                  NPCascadeTable::BachCOSPA,
                  NPCascadeTable::MassXi,
                  NPCascadeTable::MassOmega,
                  NPCascadeTable::MassV0,
                  NPCascadeTable::CascRadius,
                  NPCascadeTable::V0Radius,
                  NPCascadeTable::CascLenght,
                  NPCascadeTable::V0Lenght,
                  NPCascadeTable::CascNClusITS,
                  NPCascadeTable::CascNClusTPC,
                  NPCascadeTable::ProtonNClusTPC,
                  NPCascadeTable::PionNClusTPC,
                  NPCascadeTable::BachKaonNClusTPC,
                  NPCascadeTable::BachPionNClusTPC,
                  NPCascadeTable::ProtonTPCSignal,
                  NPCascadeTable::PionTPCSignal,
                  NPCascadeTable::BachKaonTPCSignal,
                  NPCascadeTable::BachPionTPCSignal,
                  NPCascadeTable::ProtonTPCNSigma,
                  NPCascadeTable::PionTPCNSigma,
                  NPCascadeTable::BachKaonTPCNSigma,
                  NPCascadeTable::BachPionTPCNSigma,
                  NPCascadeTable::ProtonHasTOF,
                  NPCascadeTable::PionHasTOF,
                  NPCascadeTable::BachKaonHasTOF,
                  NPCascadeTable::BachPionHasTOF,
                  NPCascadeTable::ProtonTOFSignal,
                  NPCascadeTable::PionTOFSignal,
                  NPCascadeTable::BachKaonTOFSignal,
                  NPCascadeTable::BachPionTOFSignal,
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















// // Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// // See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// // All rights not expressly granted are reserved.
// //
// // This software is distributed under the terms of the GNU General Public
// // License v3 (GPL Version 3), copied verbatim in the file "COPYING".
// //
// // In applying this license CERN does not waive the privileges and immunities
// // granted to it by virtue of its status as an Intergovernmental Organization
// // or submit itself to any jurisdiction.

// ///
// /// \file LFNonPromptCascadeTable.h
// /// \brief Non prompt cascade tables
// ///

// #include "Framework/AnalysisDataModel.h"
// #include "Framework/AnalysisTask.h"
// #include "Framework/ASoAHelpers.h"

// #ifndef PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_
// #define PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_

// namespace o2::aod
// {
// namespace NPCascadeTable /// Per cosa sta NS? // come si aggiungono nuove variabili? (c'è una qualche lista?)
// {
// DECLARE_SOA_COLUMN(Pt, pt, float);
// DECLARE_SOA_COLUMN(Eta, eta, float);
// DECLARE_SOA_COLUMN(Phi, phi, float);

// DECLARE_SOA_COLUMN(CascDCAxy, cascDCAxy, float);
// DECLARE_SOA_COLUMN(CascDCAz, cascDCAz, float);
// DECLARE_SOA_COLUMN(PtrackDCAxy, ptrackDCAxy, float);
// DECLARE_SOA_COLUMN(PtrackDCAz, ptrackDCAz, float);
// DECLARE_SOA_COLUMN(NtrackDCAxy, ntrackDCAxy, float);
// DECLARE_SOA_COLUMN(NtrackDCAz, ntrackDCAz, float);
// DECLARE_SOA_COLUMN(BachDCAxy, bachDCAxy, float);
// DECLARE_SOA_COLUMN(BachDCAz, bachDCAz, float);

// DECLARE_SOA_COLUMN(CascCosPA, casccosPA, int);
// DECLARE_SOA_COLUMN(V0COSPA, v0cosPA, int);
// DECLARE_SOA_COLUMN(BachCOSPA, bachcosPA, int);

// DECLARE_SOA_COLUMN(MassXi, massXi, float);
// DECLARE_SOA_COLUMN(MassOmega, massOmega, float);
// DECLARE_SOA_COLUMN(MassV0, massV0, float);

// DECLARE_SOA_COLUMN(CascRadius, cascRadius, float);
// DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);

// DECLARE_SOA_COLUMN(ITSnCls, itsNCls, int);
// DECLARE_SOA_COLUMN(PtrackTPCSignal, ptrackTPCSignal, float); 
// DECLARE_SOA_COLUMN(NtrackTPCSignal, ntrackTPCSignal, float); 
// DECLARE_SOA_COLUMN(BachTPCSignal, bachTPCSignal, float); 
// DECLARE_SOA_COLUMN(PtrackTPCNSigma, ptrackTPCNSigma, float); 
// DECLARE_SOA_COLUMN(NtrackTPCNSigma, ntrackTPCNSigma, float); 
// DECLARE_SOA_COLUMN(BachTPCNSigma, bachTPCNSigma, float); 
// DECLARE_SOA_COLUMN(PtrackHasTOF, ptrackhasTOF, bool); 
// DECLARE_SOA_COLUMN(NtrackHasTOF, ntrackhasTOF, bool); 
// DECLARE_SOA_COLUMN(NClusTPC, nclusTPC, int); 
// } // namespace NPCascadeTable

// DECLARE_SOA_TABLE(NonPromptCascadeTable, "AOD", "NPCASCTABLE",
//                   NPCascadeTable::Pt,
//                   NPCascadeTable::Eta,
//                   NPCascadeTable::Phi,
//                   NPCascadeTable::CascDCAxy,
//                   NPCascadeTable::CascDCAz,
//                   NPCascadeTable::PtrackDCAxy,
//                   NPCascadeTable::PtrackDCAz,
//                   NPCascadeTable::NtrackDCAxy,
//                   NPCascadeTable::NtrackDCAz,
//                   NPCascadeTable::BachDCAxy,
//                   NPCascadeTable::BachDCAz,
//                   NPCascadeTable::CascCosPA,
//                   NPCascadeTable::V0COSPA,
//                   NPCascadeTable::BachCOSPA,
//                   NPCascadeTable::MassXi,
//                   NPCascadeTable::MassOmega,
//                   NPCascadeTable::MassV0,
//                   NPCascadeTable::CascRadius,
//                   NPCascadeTable::V0Radius,
//                   NPCascadeTable::ITSnCls,
//                   NPCascadeTable::PtrackTPCSignal,
//                   NPCascadeTable::NtrackTPCSignal,
//                   NPCascadeTable::BachTPCSignal,
//                   NPCascadeTable::PtrackTPCNSigma,
//                   NPCascadeTable::NtrackTPCNSigma,
//                   NPCascadeTable::BachTPCNSigma,
//                   NPCascadeTable::PtrackHasTOF,
//                   NPCascadeTable::NtrackHasTOF,
//                   NPCascadeTable::NClusTPC)


// DECLARE_SOA_TABLE(NonPromptCascadeTableMC, "AOD", "NPCASCTABLEMC",
//                   NPCascadeTable::Pt,
//                   NPCascadeTable::Eta,
//                   NPCascadeTable::Phi,
//                   NPCascadeTable::TPCInnerParam,
//                   NPCascadeTable::Beta,
//                   NPCascadeTable::Zvertex,
//                   NPCascadeTable::DCAxy,
//                   NPCascadeTable::DCAz,
//                   NPCascadeTable::TPCsignal,
//                   NPCascadeTable::ITSchi2,
//                   NPCascadeTable::TPCchi2,
//                   NPCascadeTable::Flags,
//                   NPCascadeTable::TPCfindableCls,
//                   NPCascadeTable::TPCcrossedRows,
//                   NPCascadeTable::ITSclsMap,
//                   NPCascadeTable::TPCnCls,
//                   NPCascadeTable::gPt,
//                   NPCascadeTable::gEta,
//                   NPCascadeTable::gPhi,
//                   NPCascadeTable::PDGcode)

// } // namespace o2::aod

// #endif // PWGLF_DATAMODEL_LFNONPROMPTCASCADETABLES_H_
