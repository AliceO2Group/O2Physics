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
//
/// \file strangenessMasks.h
/// \brief Defines selection criteria and bit masks for identifying v0 and cascade particles
/// \author Roman Nepeivoda (roman.nepeivoda@cern.ch)

#ifndef PWGLF_UTILS_STRANGENESSMASKS_H_
#define PWGLF_UTILS_STRANGENESSMASKS_H_

enum SelectionsCombined : int { selV0CosPA = 0,
                                selV0Radius,
                                selV0RadiusMax,
                                selDCANegToPV,
                                selDCAPosToPV,
                                selDCAV0Dau,
                                selK0ShortRapidity,
                                selLambdaRapidity,
                                selTPCPIDPositivePion,
                                selTPCPIDNegativePion,
                                selTPCPIDPositiveProton,
                                selTPCPIDNegativeProton,
                                selTOFDeltaTPositiveProtonLambda,
                                selTOFDeltaTPositivePionLambda,
                                selTOFDeltaTPositivePionK0Short,
                                selTOFDeltaTNegativeProtonLambda,
                                selTOFDeltaTNegativePionLambda,
                                selTOFDeltaTNegativePionK0Short,
                                selTOFNSigmaPositiveProtonLambda, // Nsigma
                                selTOFNSigmaPositivePionLambda,   // Nsigma
                                selTOFNSigmaPositivePionK0Short,  // Nsigma
                                selTOFNSigmaNegativeProtonLambda, // Nsigma
                                selTOFNSigmaNegativePionLambda,   // Nsigma
                                selTOFNSigmaNegativePionK0Short,  // Nsigma
                                selK0ShortCTau,
                                selLambdaCTau,
                                selK0ShortArmenteros,
                                selPosGoodTPCTrack, // at least min # TPC rows
                                selNegGoodTPCTrack, // at least min # TPC rows
                                selPosGoodITSTrack, // at least min # ITS clusters
                                selNegGoodITSTrack, // at least min # ITS clusters
                                selPosItsOnly,
                                selNegItsOnly,
                                selPosNotTPCOnly,
                                selNegNotTPCOnly,
                                selConsiderK0Short,    // for mc tagging
                                selConsiderLambda,     // for mc tagging
                                selConsiderAntiLambda, // for mc tagging
                                selPhysPrimK0Short,    // for mc tagging
                                selPhysPrimLambda,     // for mc tagging
                                selPhysPrimAntiLambda, // for mc tagging
                                selPosEta,
                                selNegEta,
                                selDauDCA,
                                // cascade selections
                                selCascCosPA,
                                selMassWinXi,
                                selMassWinOmega,
                                selDCACascDau,
                                selDCAV0ToPV,
                                selCascRadius,
                                selCascRadiusMax,
                                selBachBaryon,
                                selRejCompXi,
                                selRejCompOmega,
                                selBachToPV,
                                selMesonToPV,
                                selBaryonToPV,
                                selXiRapidity,
                                selXiCTau,
                                selConsiderXi,
                                selConsiderAntiXi,
                                selPhysPrimXi,
                                selPhysPrimAntiXi,
                                selOmegaRapidity,
                                selOmegaCTau,
                                selConsiderOmega,
                                selConsiderAntiOmega,
                                selPhysPrimOmega,
                                selPhysPrimAntiOmega,
                                selBachItsOnly,
                                selBachNotTPCOnly,
                                selBachGoodTPCTrack,
                                selBachGoodITSTrack,
                                selTPCPIDBachPion,
                                selTPCPIDBachKaon,
                                selTOFNSigmaBachPionXi,
                                selTOFNSigmaBachKaonOmega,
                                selTOFNSigmaPositiveProtonLambdaXi,
                                selTOFNSigmaPositiveProtonLambdaOmega,
                                selTOFNSigmaPositivePionLambdaXi,
                                selTOFNSigmaPositivePionLambdaOmega,
                                selTOFNSigmaNegativePionLambdaXi,
                                selTOFNSigmaNegativePionLambdaOmega,
                                selTOFNSigmaNegativeProtonLambdaXi,
                                selTOFNSigmaNegativeProtonLambdaOmega,
                                selTOFDeltaTPositiveProtonLambdaXi,
                                selTOFDeltaTPositiveProtonLambdaOmega,
                                selTOFDeltaTPositivePionLambdaXi,
                                selTOFDeltaTPositivePionLambdaOmega,
                                selTOFDeltaTNegativeProtonLambdaXi,
                                selTOFDeltaTNegativeProtonLambdaOmega,
                                selTOFDeltaTNegativePionLambdaXi,
                                selTOFDeltaTNegativePionLambdaOmega,
                                selTOFDeltaTBachPionXi,
                                selTOFDeltaTBachKaonOmega,
                                selBachEta,
                                selLambdaMassWin,
                                selCount,
};

static constexpr int kSelNum = static_cast<int>(SelectionsCombined::selCount);

// constants
const float ctauxiPDG = 4.91;     // from PDG
const float ctauomegaPDG = 2.461; // from PDG

// bit masks
std::bitset<kSelNum> maskTopologicalV0;
std::bitset<kSelNum> maskTopologicalCasc;

std::bitset<kSelNum> maskKinematicV0;
std::bitset<kSelNum> maskKinematicCasc;

std::bitset<kSelNum> maskTrackPropertiesV0;
std::bitset<kSelNum> maskTrackPropertiesCasc;

std::bitset<kSelNum> maskK0ShortSpecific;
std::bitset<kSelNum> maskLambdaSpecific;
std::bitset<kSelNum> maskAntiLambdaSpecific;
std::bitset<kSelNum> maskXiSpecific;
std::bitset<kSelNum> maskAntiXiSpecific;
std::bitset<kSelNum> maskOmegaSpecific;
std::bitset<kSelNum> maskAntiOmegaSpecific;

std::bitset<kSelNum> maskSelectionK0Short;
std::bitset<kSelNum> maskSelectionLambda;
std::bitset<kSelNum> maskSelectionAntiLambda;
std::bitset<kSelNum> maskSelectionXi;
std::bitset<kSelNum> maskSelectionAntiXi;
std::bitset<kSelNum> maskSelectionOmega;
std::bitset<kSelNum> maskSelectionAntiOmega;

#endif // PWGLF_UTILS_STRANGENESSMASKS_H_
