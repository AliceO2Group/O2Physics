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

/// \file CorrelatedBkgTagger.h
/// \brief Definitions of channels for correlated bkgs
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic of Turin and INFN

#ifndef PWGHF_CORE_CORRELATEDBKGTAGGER_H_
#define PWGHF_CORE_CORRELATEDBKGTAGGER_H_

namespace o2::hf_corrbkg
{

    // D± → K± K∓ π±
    enum DecayChannelResoDplus {
        PhiPi = 0,
        K0starK,
        // K01430K,
        // RhoPi,
        // f2_1270Pi
    };

    // Ds± → K± K∓ π±
    enum DecayChannelResoDs {
        PhiPi = DecayChannelResoDplus::PhiPi,
        K0starK,
        // PhiRho,
        // f2_1270Pi,
        // K0starPi,
        // RhoK,
        // EtaPi
    };

    // enum DecayChannelResoD0 {
    //     PhiPi = DecayChannelResoDs::PhiPi,
    //     K0starPi,
    //     RhoK,
    //     K0starPi0,
    //     K0starPi,
    //     RhoPi,
    // };

    // enum DecayChannelResoLc {
    //     K0starP = DecayChannelResoD0::PhiPi,
    //     DeltaK,
    //     Lambda1520K,
    //     PhiP,
    // };

    // enum DecayChannelResoXic {
    //     K0starP = DecayChannelResoLc::K0starP,
    //     PhiP,
    // };

    // Function to associate vectors with enum values
    std::vector<int> getFinalDaughters(DecayType finalState) {
        switch(finalState) {
            case DecayChannelCorrBkg::DstarToPiKPiBkg:
                return std::array{+kPiPlus, +kPiPlus, -kKPlus};
            default:
                return std::array{};
        }
    }
} // namespace o2::hf_corrbkg

#endif // PWGHF_CORE_CORRELATEDBKGTAGGER_H_
