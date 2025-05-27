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

/// \file CorrelatedBkgs.h
/// \brief Definitions of channels for correlated bkgs
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic of Turin and INFN
/// \author Stefano Politanò <stefano.politano@cern.ch>, CERN

#include <unordered_map>
#include <array>
#include <variant>

#include "/home/mdicosta/alice/O2Physics/Common/Core/RecoDecay.h"

#ifndef PWGHF_CORE_CORRELATEDBKGS_H_
#define PWGHF_CORE_CORRELATEDBKGS_H_

using namespace o2::constants::physics;

namespace o2::hf_corrbkg
{

    enum FinalStates2Prongs {
        PiK = 1,
        KPi,
        KPiPi0,
        PiPi,
        PiPiPi0,
    };

    std::unordered_map<FinalStates2Prongs, std::vector<int> > finalStates2Prongs = 
    {
        {FinalStates2Prongs::PiK, std::vector<int>{+kPiPlus, -kKPlus}},
        {FinalStates2Prongs::KPi, std::vector<int>{-kKPlus, +kPiPlus}},
        {FinalStates2Prongs::KPiPi0, std::vector<int>{-kKPlus, +kPiPlus, +kPi0}},
        {FinalStates2Prongs::PiPi, std::vector<int>{+kPiMinus, +kPiPlus}},
        {FinalStates2Prongs::PiPiPi0, std::vector<int>{+kPiMinus, +kPiPlus, +kPi0}}
    };

    enum FinalStates3Prongs {
        KKPi = 1,
        KKPiPi0,
        KPiPi,
        KPiPiPi0,
        PiPiPi,
        PiPiPiPi0,
        ProtonPiPi,
        ProtonKPi,
        NFinalStates
    };

    std::unordered_map<FinalStates3Prongs, std::vector<int> > finalStates3Prongs = 
    {
        {FinalStates3Prongs::KKPi, std::vector<int>{-kKPlus, +kKPlus, +kPiPlus}},
        {FinalStates3Prongs::KKPiPi0, std::vector<int>{-kKPlus, +kKPlus, +kPiPlus, +kPi0}},
        {FinalStates3Prongs::KPiPi, std::vector<int>{+kPiPlus, -kKPlus, +kPiPlus}},
        {FinalStates3Prongs::KPiPiPi0, std::vector<int>{+kPiPlus, -kKPlus, +kPiPlus, +kPi0}},
        {FinalStates3Prongs::PiPiPi, std::vector<int>{+kPiMinus, +kPiPlus, +kPiPlus}},
        {FinalStates3Prongs::PiPiPiPi0, std::vector<int>{+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
        {FinalStates3Prongs::ProtonPiPi, std::vector<int>{kProton, -kPiPlus, +kPiPlus}},
        {FinalStates3Prongs::ProtonKPi, std::vector<int>{kProton, -kKPlus, +kPiPlus}}
    };

    // Dstar → K± K∓ π±
    namespace D0 {
      enum DecayChannelResoD0 {
        RhoK = 1,
        K0starPi0,
        K0starPiPlus,
        RhoPi,
      };

      std::unordered_map<DecayChannelResoD0, std::array<int, 2> > resoStatesD0 = 
      {
        {DecayChannelResoD0::RhoK, std::array<int, 2>{213, -kKPlus}},
        {DecayChannelResoD0::K0starPi0, std::array<int, 2>{-kK0Star892, -kKPlus}},
        {DecayChannelResoD0::K0starPiPlus, std::array<int, 2>{-323, kKPlus}},
        {DecayChannelResoD0::RhoPi, std::array<int, 2>{213, -kPiPlus}},
      };
    }

    // D± → K± K∓ π±
    namespace DPlus {
      enum DecayChannelResoDplus {
        K0starK = 1,
        K_1430K,
        PhiPi,
        RhoPi,
        f2_1270Pi,
      };

      std::unordered_map<DecayChannelResoDplus, std::array<int, 2> > resoStatesDPlus = 
      {
        {DecayChannelResoDplus::K0starK, std::array<int, 2>{-kK0Star892, -kKPlus}},
        {DecayChannelResoDplus::K_1430K, std::array<int, 2>{-10321, -kKPlus}},
        {DecayChannelResoDplus::PhiPi, std::array<int, 2>{+kPhi, +kPiPlus}},
        {DecayChannelResoDplus::RhoPi, std::array<int, 2>{+213, +kPiPlus}},
        {DecayChannelResoDplus::f2_1270Pi, std::array<int, 2>{225, +kPiPlus}},
      };
    }
    
    // Ds± → K± K∓ π±
    namespace DS {
      enum DecayChannelResoDs {
          PhiPi = 1,
          K0starK,
          PhiRho,
          f2_1270Pi,
          K0starPi,
          RhoK,
          EtaPi,
      };
        
      std::unordered_map<DecayChannelResoDs, std::array<int, 2> > resoStatesDs = 
      {
        {DecayChannelResoDs::PhiPi, std::array<int, 2>{+kPhi, +kPiPlus}},
        {DecayChannelResoDs::K0starK, std::array<int, 2>{-kK0Star892, -kKPlus}},
        {DecayChannelResoDs::PhiRho, std::array<int, 2>{-kPhi, -113}},
        {DecayChannelResoDs::f2_1270Pi, std::array<int, 2>{225, +kPiPlus}},
        {DecayChannelResoDs::K0starPi, std::array<int, 2>{-kK0Star892, -kPiPlus}},
        {DecayChannelResoDs::RhoK, std::array<int, 2>{-113, -kKPlus}},
        {DecayChannelResoDs::EtaPi, std::array<int, 2>{-221, -kPiPlus}},
      };
    }

    // Dstar → K± K∓ π±
    namespace DStar {
      enum DecayChannelResoDStarD0 {
        RhoK = 1,
        K0starPi0,
        K0starPiPlus,
        RhoPi,
      };

      std::unordered_map<DecayChannelResoDStarD0, std::array<int, 2> > resoStatesDStarD0 = 
      {
        {DecayChannelResoDStarD0::RhoK, std::array<int, 2>{213, -kKPlus}},
        {DecayChannelResoDStarD0::K0starPi0, std::array<int, 2>{-kK0Star892, -kKPlus}},
        {DecayChannelResoDStarD0::K0starPiPlus, std::array<int, 2>{-323, kKPlus}},
        {DecayChannelResoDStarD0::RhoPi, std::array<int, 2>{213, -kPiPlus}},
      };
    }

    // Lc → p K∓ π±
    namespace LambdaC {
      enum DecayChannelResoLambdaC {
        K0starP = 1,
        DeltaK,
        Lambda1520Pi,
        ProtonPhi,
      }; 

      std::unordered_map<DecayChannelResoLambdaC, std::array<int, 2> > resoStatesLambdaC = 
      {
        {DecayChannelResoLambdaC::K0starP, std::array<int, 2>{-kK0Star892, -kProton}},
        {DecayChannelResoLambdaC::DeltaK, std::array<int, 2>{+2224, -kKPlus}},
        {DecayChannelResoLambdaC::Lambda1520Pi, std::array<int, 2>{+102134, -kPiPlus}},
        {DecayChannelResoLambdaC::ProtonPhi, std::array<int, 2>{+kProton, -kPhi}},
      };
    }
    
    // Xic → p K∓ π±
    namespace XiC {
      enum DecayChannelResoXiC {
        K0starP = 1,
        ProtonPhi,
        SigmaPiPi,
      };
      
      std::unordered_map<DecayChannelResoXiC, std::array<int, 2> > resoStatesXiC = 
      {
        {DecayChannelResoXiC::K0starP, std::array<int, 2>{-kK0Star892, -kProton}},
        {DecayChannelResoXiC::ProtonPhi, std::array<int, 2>{+kProton, -kPhi}},
        {DecayChannelResoXiC::SigmaPiPi, std::array<int, 2>{+kProton, -kPhi}},
      };
    }

    bool checkResonantDecay(std::vector<int> arrDaughIndex, std::array<int, 2> arrPDGResonant) {
      std::cout << "###";
      for (int iArrDaugh : arrDaughIndex) {
        bool findDaug = false;
        for (int i = 0; i < 2; ++i) {
          std::cout << "Checking daughter PDG: " << iArrDaugh << " against resonant PDG: " << arrPDGResonant[i] << std::endl;
          if (std::abs(iArrDaugh) == arrPDGResonant[i]) {
            findDaug = true;
          }
        }
        if (!findDaug) {
          std::cout << "###";
          return false; // If any daughter does not match, return false
        }
      }
      std::cout << "###";
      return true;
    }

    /// Check if the decay is resonant
    /// \tparam arrDaughIndex index of the particle daughters at resonance level
    /// \tparam arrPDGResonant PDG code of the resonant decay
    /// \return true if the decay is resonant
    void flagResonantDecay(int motherPdg, int8_t* channel, std::vector<int> arrDaughIndex) {
      switch (motherPdg) {
        case Pdg::kD0:
          for (const auto& [flag, pdgCodes] : D0::resoStatesD0) {
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = (1 << flag);
              break;
            }
          }
          break;
        case Pdg::kDPlus:
          for (const auto& [flag, pdgCodes] : DPlus::resoStatesDPlus) {
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = (1 << flag);
              break;
            }
          }
          break;
        case Pdg::kDS:
          for (const auto& [flag, pdgCodes] : DS::resoStatesDs) {
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = (1 << flag);
              break;
            }
          }
          break;
        case Pdg::kDStar:
          for (const auto& [flag, pdgCodes] : DStar::resoStatesDStarD0) {
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = (1 << flag);
              break;
            }
          }
          break;
        case Pdg::kLambdaCPlus:
          for (const auto& [flag, pdgCodes] : LambdaC::resoStatesLambdaC) {
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = (1 << flag);
              break;
            }
          }
          break;
        case Pdg::kXiCPlus:
          for (const auto& [flag, pdgCodes] : XiC::resoStatesXiC) {
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = (1 << flag);
              break;
            }
          }
          break;
      }
    }

    template <bool matchKinkedDecayTopology = false, bool matchInteractionsWithMaterial = false>
    int matchFinalStateCorrBkgs(int pdgMother, aod::McParticles const& mcParticles, auto arrayDaughters, int8_t* flag, int8_t* sign, int8_t* channel, int depth, int8_t* nKinkedTracks, int8_t* nInteractionsWithMaterial) {
      std::cout << "Entering flag value: " << static_cast<int>(*flag) << std::endl;
      std::vector<int> arrResoDaughIndex = {};
      int indexRec = -1; // Index of the matched reconstructed candidate
      // for (std::size_t iProng = 0; iProng < 3; ++iProng) {
      //   if (!arrayDaughters[iProng].has_mcParticle()) {
      //       return -1;
      //   }
      //   auto particleI = arrayDaughters[iProng].mcParticle(); // ith daughter particle
      //   auto motherI = particleI.template mothers_first_as<aod::McParticles>();
      //   auto pdgI = particleI.pdgCode();
      //   auto pdgMotherI = motherI.pdgCode();
      //   auto pdgMotherII = -1;
      //   LOG(info) << "Daughter " << iProng << " PDG: " << pdgI << " motherI: " << pdgMotherI;
      // }
      LOG(info) << "Matching daughters ... ";
      int matchedMotherId = 0;
      
      // auto finalStates = finalStates3Prongs;
      // if (std::abs(pdgMother) == Pdg::kD0) {
      //   finalStates = finalStates2Prongs;
      // }

      for (const auto& [chn, finalState] : finalStates3Prongs) {
        if (finalState.size() == 4) {                     // Partly Reco 3-prong decays
          std::array<int, 4> finalStateParts = std::array{finalState[0], finalState[1], finalState[2], finalState[3]};
          if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
          } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks);
          } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nullptr, nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, false, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth);
          }
        } else if (finalState.size() == 3 && std::abs(pdgMother) != Pdg::kD0) { // Fully Reco 3-prong decays
          std::array<int, 3> finalStateParts = std::array{finalState[0], finalState[1], finalState[2]};
          if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
          } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks);
          } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nullptr, nInteractionsWithMaterial);
          } else {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, false>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth);
          }
        } else if (finalState.size() == 3 && std::abs(pdgMother) == Pdg::kD0) {  // Partly Reco 2-prong decays
          std::array<int, 2> finalStateParts = std::array{finalState[0], finalState[1]};
          if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
          } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks);
          } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nullptr, nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign);
          }
        } else if (finalState.size() == 2) { // Fully Reco 2-prong decays
          std::array<int, 2> finalStateParts = std::array{finalState[0], finalState[1]};
          if constexpr (matchKinkedDecayTopology && matchInteractionsWithMaterial) { 
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks, nInteractionsWithMaterial);
          } else if constexpr (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nKinkedTracks);
          } else if constexpr (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign, depth, nullptr, nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false>(mcParticles, arrayDaughters, pdgMother, finalStateParts, true, sign);
          }
        } else {
          LOG(info) << "Final state size not supported: " << finalState.size();
          continue; // Skip unsupported final states
        }
        LOG(info) << "IndexRec: " << indexRec;
        if (indexRec > -1) {
          std::cout << "Matched final state: " << chn << " with PDG code: " << pdgMother << std::endl;
          switch (pdgMother) {
            case Pdg::kD0:
              *flag = (*sign) * chn + 10;
              break;
            case Pdg::kDPlus:
              *flag = (*sign) * chn + 20;
              break;
              case Pdg::kDS:
              *flag = (*sign) * chn + 30;
              break;
            case Pdg::kDStar:
              *flag = (*sign) * chn + 40;
              break;
            case Pdg::kLambdaCPlus:
              *flag = (*sign) * chn + 50;
              break;
            case Pdg::kXiCPlus:
              *flag = (*sign) * chn + 60;
              break;
            default:
              LOG(info) << "Unknown mother PDG code: " << pdgMother << ", skipping.";
              continue; // Skip unknown mother PDG codes
          }
          
          // Flag the resonant decay channel
          int resoMaxDepth = 1;
          int NDaughtersResonant = 2;
          if (std::abs(pdgMother) == Pdg::kDStar) { 
            resoMaxDepth = 2; // Flag D0 resonances 
          }
          RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrResoDaughIndex, std::array{0}, resoMaxDepth);
          std::vector<int> arrPDGDaugh = {};
          if (arrResoDaughIndex.size() == NDaughtersResonant) {
            for (auto iProng = 0u; iProng < arrResoDaughIndex.size(); ++iProng) {
              auto daughI = mcParticles.rawIteratorAt(arrResoDaughIndex[iProng]);
              if ( (std::abs(pdgMother) == Pdg::kDStar || std::abs(pdgMother) == Pdg::kXiCPlus) 
                    && std::abs(daughI.pdgCode() == kPiPlus && arrPDGDaugh.size() >= 2)) {
                continue; // Skip the pion from D* decay and the second pion from XiC --> Sigma Pi Pi
              }
              arrPDGDaugh.push_back(std::abs(daughI.pdgCode()));
            }
            flagResonantDecay(pdgMother, channel, arrPDGDaugh);
            LOG(info) << "[matchFinalStateDMeson] Flag set to: " << static_cast<int>(*flag) << " sign: " << static_cast<int>(*sign) << " for channel: " << chn;
            break;
          }
        }
      }
      return indexRec;
    }

} // namespace o2::hf_corrbkg

#endif // PWGHF_CORE_CORRELATEDBKGS_H_
