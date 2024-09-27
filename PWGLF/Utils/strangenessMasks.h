// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
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
/// \author Roman Nepeivoda (roman.nepeivoda@cern.ch)
/// \since August 27, 2024

#ifndef PWGLF_UTILS_STRANGENESSMASKS_H_
#define PWGLF_UTILS_STRANGENESSMASKS_H_

enum selectionsCombined : int { selV0CosPA = 0,
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

static constexpr int selNum = static_cast<int>(selectionsCombined::selCount);

// constants
const float ctauxiPDG = 4.91;     // from PDG
const float ctauomegaPDG = 2.461; // from PDG

// bit masks
std::bitset<selNum> maskTopologicalV0;
std::bitset<selNum> maskTopologicalCasc;

std::bitset<selNum> maskKinematicV0;
std::bitset<selNum> maskKinematicCasc;

std::bitset<selNum> maskTrackPropertiesV0;
std::bitset<selNum> maskTrackPropertiesCasc;

std::bitset<selNum> maskK0ShortSpecific;
std::bitset<selNum> maskLambdaSpecific;
std::bitset<selNum> maskAntiLambdaSpecific;
std::bitset<selNum> maskXiSpecific;
std::bitset<selNum> maskAntiXiSpecific;
std::bitset<selNum> maskOmegaSpecific;
std::bitset<selNum> maskAntiOmegaSpecific;

std::bitset<selNum> maskSelectionK0Short;
std::bitset<selNum> maskSelectionLambda;
std::bitset<selNum> maskSelectionAntiLambda;
std::bitset<selNum> maskSelectionXi;
std::bitset<selNum> maskSelectionAntiXi;
std::bitset<selNum> maskSelectionOmega;
std::bitset<selNum> maskSelectionAntiOmega;

#endif // PWGLF_UTILS_STRANGENESSMASKS_H_
