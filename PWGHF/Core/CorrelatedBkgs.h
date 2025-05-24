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

#include <unordered_map>
#include <array>
#include <variant>

#ifndef PWGHF_CORE_CORRELATEDBKGTAGGER_H_
#define PWGHF_CORE_CORRELATEDBKGTAGGER_H_

using namespace o2::constants::physics;

namespace o2::hf_corrbkg
{

    enum FinalStatesDMesons {
        KKPi = 0,
        KPiPi,
        KPiPiPi0,
        PiPiPi,
        PiPiPiPi0,
        NFinalStatesDMesons
    };

    // D± → K± K∓ π±
    enum DecayChannelResoDplus {
        PhiPiDplus = 0,
        K0starKDplus,
        // K01430K,
        // RhoPi,
        // f2_1270Pi
    };

    // Ds± → K± K∓ π±
    enum DecayChannelResoDs {
        PhiPiDs = 0,
        K0starKDs,
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


    // using FinalStateVariant = std::variant<
    //     std::array<int, 3>,
    //     std::array<int, 4>
    // >;

    // std::unordered_map<int, FinalStateVariant> finalStates = 
    // {
    //     {FinalStatesDMesons::KKPi, std::array<int,3>{+kKPlus, +kKPlus, +kPiPlus}},
    //     {FinalStatesDMesons::KPiPiPi0, std::array<int,4>{+kKPlus, +kKPlus, +kPiPlus, +kPi0}}
    // };

    std::unordered_map<FinalStatesDMesons, std::pair<std::array<int,3>, bool>> finalStates = 
    {
        {FinalStatesDMesons::KKPi, {std::array<int,3>{-kKPlus, +kKPlus, +kPiPlus}, false}},
        {FinalStatesDMesons::KPiPi, {std::array<int,3>{+kPiPlus, -kKPlus, +kPiPlus}, false}},
        {FinalStatesDMesons::KPiPiPi0, {std::array<int,3>{+kPiPlus, -kKPlus, +kPiPlus}, true}},
        {FinalStatesDMesons::PiPiPi, {std::array<int,3>{+kPiMinus, +kPiPlus, +kPiPlus}, false}},
        {FinalStatesDMesons::PiPiPiPi0, {std::array<int,3>{+kPiMinus, +kPiPlus, +kPiPlus}, true}},
    };


    template <bool matchKinkedDecayTopology = false, bool matchInteractionsWithMaterial = false>
    int matchFinalStateDMeson(aod::McParticles const& mcParticles, auto arrayDaughters, int pdg, int8_t* flag, int8_t* sign, int depth, int8_t* nKinkedTracks, int8_t* nInteractionsWithMaterial) {
        int indexRec = -1;
        LOG(info) << "matchFinalStateDMeson: " << pdg;
        for (const std::pair<int, std::pair<std::array<int,3>, bool>>& finalState : finalStates) {
            if (finalState.second.second) { // 4 prong decay, partly reconstructed
                LOG(info) << "[4prong] Taking channel with final state: " << finalState.second.first[0] << " " << finalState.second.first[1] << " " << finalState.second.first[2] << " " << finalState.second.first[3] << ", has 4 daughters: " << finalState.second.second;
                LOG(info) << "[4prong] Trying 4-prong decay ... ";
                LOG(info) << "[4prong] Matching channel with final state: " << finalState.second.first[0] << " " << finalState.second.first[1] << " " << finalState.second.first[2];
                if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, true, true, true>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
                } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, true, true, false>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, sign, depth, nKinkedTracks);
                } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, true, false, true>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, sign, depth, nullptr, nInteractionsWithMaterial);
                } else {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, true, false, false>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, &sign, depth);
                }
            } else {
                LOG(info) << "[3prong] Trying 3-prong decay ... ";
                LOG(info) << "[3prong] Taking channel with final state: " << finalState.second.first[0] << " " << finalState.second.first[1] << " " << finalState.second.first[2] << ", has 4 daughters: " << finalState.second.second;
                if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
                } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, sign, depth, nKinkedTracks);
                } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, sign, depth, nullptr, nInteractionsWithMaterial);
                } else {
                    indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, false>(mcParticles, arrayDaughters, pdg, finalState.second.first, true, &sign, depth);
                }
            }
            if (indexRec > -1) {
                LOG(info) << "Found channel with final state: " << finalState.second.first[0] << " " << finalState.second.first[1] << " " << finalState.second.first[2] << ", has 4 prongs: " << finalState.second.second;
                *flag = (*sign) * (1 << finalState.first);
                break;
            }
        }
        return indexRec;
    }





    // Function to associate vectors with enum values
    // std::vector<int> getFinalDaughters(DecayType finalState) {
    //     switch(finalState) {
    //         case DecayChannelCorrBkg::DstarToPiKPiBkg:
    //             return std::array{+kPiPlus, +kPiPlus, -kKPlus};
    //         default:
    //             return std::array{};
    //     }
    // }

    // std::vector<Pdg> getFinalStateDMeson(int index) {
    //     switch (index) {
    //         case KKPi:
    //             return std::vector<Pdg>{+kKPlus, -kKPlus, +kPiPlus};
    //         case KPiPi:
    //             return std::vector<Pdg>{+kKPlus, +kPiPlus, +kPiPlus};
    //         case KPiPiPi0:
    //             return std::vector<Pdg>{+kKPlus, +kPiPlus, +kPiPlus, +kPi0};
    //         case PiPiPi:
    //             return std::vector<Pdg>{+kPiPlus, +kPiPlus, +kPiPlus};
    //         case PiPiPiPi0:
    //             return std::vector<Pdg>{+kPiPlus, +kPiPlus, +kPiPlus, +kPi0};
    //         default:
    //             throw std::out_of_range("Invalid index for final state D meson");
    //     }
    // }


    // template <int N>
    // std::array<int, N> getFinalStateDMeson(int index) {
    //     if constexpr (N == 3) {
    //         switch (index) {
    //             case FinalStatesDMesons::KKPi:
    //                 return std::array<int, N>{+kKPlus, -kKPlus, +kPiPlus};
    //             case FinalStatesDMesons::KPiPi:
    //                 return std::array<int, N>{+kKPlus, +kPiPlus, +kPiPlus};
    //             case FinalStatesDMesons::PiPiPi:
    //                 return std::array<int, N>{+kPiPlus, +kPiPlus, +kPiPlus};
    //             default:
    //                 throw std::out_of_range("Invalid index for final state D meson");
    //         }
    //     } else if constexpr (N == 4) {
    //         switch (index) {
    //             case FinalStatesDMesons::KPiPiPi0:
    //                 return std::array<int, N>{+kKPlus, +kPiPlus, +kPiPlus, +kPi0};
    //             case FinalStatesDMesons::PiPiPiPi0:
    //                 return std::array<int, N>{+kPiPlus, +kPiPlus, +kPiPlus, +kPi0};
    //             default:
    //                 throw std::out_of_range("Invalid index for final state D meson");
    //         }
    //     } else {
    //         throw std::invalid_argument("Invalid size for final state D meson");
    //     }
    // }

    // int getFinalStateParts(int index) {
    //     switch (index) {
    //         case FinalStatesDMesons::KKPi:
    //             return 3;
    //         case FinalStatesDMesons::KPiPi:
    //             return 3;
    //         case FinalStatesDMesons::KPiPiPi0:
    //             return 4;
    //         case FinalStatesDMesons::PiPiPi:
    //             return 3;
    //         case FinalStatesDMesons::PiPiPiPi0:
    //             return 4;
    //         default:
    //             throw std::out_of_range("Invalid index for final state D meson");
    //     }
    // }

    // template <bool matchKinkedDecayTopology = false, bool matchInteractionsWithMaterial = false>
    // int matchFinalStateDMeson(int iChannel, aod::McParticles const& mcParticles, auto arrayDaughters, int pdg, int8_t* sign, int depth, int8_t* nKinkedTracks, int8_t* nInteractionsWithMaterial) {
    //     int indexRec = -1;
    //     int nFinalStateParts = getFinalStateParts(iChannel);
    //     if (nFinalStateParts == 3) {
    //         auto finalState = getFinalStateDMeson<3>(iChannel);
    //         if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
    //             indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, pdg, finalState, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
    //         } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
    //             indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, pdg, finalState, true, sign, depth, nKinkedTracks);
    //         } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
    //             indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, pdg, finalState, true, sign, depth, nullptr, nInteractionsWithMaterial);
    //         } else {
    //             indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdg, finalState, true, &sign, depth);
    //         }
    //         if (indexRec > -1) {
    //             LOG(info) << "Matched final state: " << iChannel;
    //             return indexRec;
    //         }
    //     } else if (nFinalStateParts == 4) {
    //         auto finalState = getFinalStateDMeson<4>(iChannel);
    //         if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
    //             indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, pdg, finalState, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
    //         } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
    //             indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, pdg, finalState, true, sign, depth, nKinkedTracks);
    //         } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
    //             indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, pdg, finalState, true, sign, depth, nullptr, nInteractionsWithMaterial);
    //         } else {
    //             indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdg, finalState, true, &sign, depth);
    //         }
    //         if (indexRec > -1) {
    //             LOG(info) << "Matched final state: " << iChannel;
    //             return indexRec;
    //         }
    //     } else {
    //         throw std::invalid_argument("Invalid size for final state D meson");
    //     }
    // }

} // namespace o2::hf_corrbkg

#endif // PWGHF_CORE_CORRELATEDBKGTAGGER_H_
