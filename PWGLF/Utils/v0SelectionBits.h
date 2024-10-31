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

#ifndef PWGLF_UTILS_V0SELECTIONBITS_H_
#define PWGLF_UTILS_V0SELECTIONBITS_H_

namespace v0data
{
// provides simple switches
enum selection : uint64_t { selCosPA = 0,
                            selRadius,
                            selRadiusMax,
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
};
} // namespace v0data

#endif // PWGLF_UTILS_V0SELECTIONBITS_H_
