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
        KPi = 1,
        KK,
        KMinusPiPi0,
        PiPi,
        PiPiPi0,
        NFinalStates2P
    };

    std::unordered_map<FinalStates2Prongs, std::vector<int> > finalStates2Prongs = 
    {
        {FinalStates2Prongs::KPi,          std::vector<int>{+kKMinus,   +kPiPlus}},
        {FinalStates2Prongs::KK,           std::vector<int>{+kKMinus,   +kKPlus}},
        {FinalStates2Prongs::KMinusPiPi0,  std::vector<int>{+kKMinus,   +kPiPlus, +kPi0}},
        {FinalStates2Prongs::PiPi,         std::vector<int>{+kPiMinus,  +kPiPlus}},
        {FinalStates2Prongs::PiPiPi0,      std::vector<int>{+kPiMinus,  +kPiPlus, +kPi0}}
    };

    enum FinalStates3Prongs {
        KKPi = 1,
        KKPiPi0,
        KMinusPiPi,
        KMinusPiPiPi0,
        KPlusPiPi,
        KPlusPiPiPi0,
        PiPiPi,
        PiPiPiPi0,
        ProtonKPi,
        ProtonKPiPi0,
        ProtonPiPi,
        ProtonKK,
        SigmaPiPi,
        NFinalStates3P
    };

    std::unordered_map<FinalStates3Prongs, std::vector<int> > finalStates3Prongs = 
    {
        {FinalStates3Prongs::KKPi,          std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus}},
        {FinalStates3Prongs::KKPiPi0,       std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus, +kPi0}},
        {FinalStates3Prongs::KMinusPiPi,    std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus}},
        {FinalStates3Prongs::KMinusPiPiPi0, std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus, +kPi0}},
        {FinalStates3Prongs::KPlusPiPi,     std::vector<int>{+kPlus,   +kPiPlus,  +kPiMinus}},
        {FinalStates3Prongs::KPlusPiPiPi0,  std::vector<int>{+kPlus,   +kPiPlus,  +kPiMinus, +kPi0}},
        {FinalStates3Prongs::PiPiPi,        std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus}},
        {FinalStates3Prongs::PiPiPiPi0,     std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus, +kPi0}},
        {FinalStates3Prongs::ProtonKPi,     std::vector<int>{+kProton,   +kKMinus,  +kPiPlus}},
        {FinalStates3Prongs::ProtonKPiPi0,  std::vector<int>{+kProton,   +kKMinus,  +kPiPlus, +kPi0}},
        {FinalStates3Prongs::ProtonPiPi,    std::vector<int>{+kProton,   +kPiMinus, +kPiPlus}},
        {FinalStates3Prongs::ProtonKK,      std::vector<int>{+kProton,   +kKMinus,  +kKPlus}},
        {FinalStates3Prongs::SigmaPiPi,     std::vector<int>{+kSigmaPlus,   +kPiMinus,  +kPiPlus}},
    };

    // Dstar → K± K∓ π±
    namespace D0 {
      enum DecayChannelResoD0 {
        RhoPi = 1,
        RhoK,
        K0starPi0,
        K0starPiPlus,
      };

      std::unordered_map<DecayChannelResoD0, std::array<int, 2> > resoStatesD0 = 
      {
        {DecayChannelResoD0::RhoPi,         std::array<int, 2>{213,            +kPiMinus}},
        {DecayChannelResoD0::RhoK,          std::array<int, 2>{213,            +kKMinus}},
        {DecayChannelResoD0::K0starPi0,     std::array<int, 2>{-kK0Star892,    +kPi0}},
        {DecayChannelResoD0::K0starPiPlus,  std::array<int, 2>{-kKPlusStar892, +kPiPlus}},
      };
    }

    // D± → K± K∓ π±
    namespace DPlus {
      enum DecayChannelResoDplus {
        PhiPi = 1,
        K0starK,
        K_1430K,
        RhoPi,
        f2_1270Pi,
      };

      std::unordered_map<DecayChannelResoDplus, std::array<int, 2> > resoStatesDPlus = 
      {
        {DecayChannelResoDplus::PhiPi,      std::array<int, 2>{+kPhi,       +kPiPlus}},
        {DecayChannelResoDplus::K0starK,    std::array<int, 2>{-kK0Star892, +kKPlus}},
        {DecayChannelResoDplus::K_1430K,    std::array<int, 2>{+10311,      +kKPlus}},
        {DecayChannelResoDplus::RhoPi,      std::array<int, 2>{+113,        +kPiPlus}},
        {DecayChannelResoDplus::f2_1270Pi,  std::array<int, 2>{+225,        +kPiPlus}},
      };

      std::unordered_map<FinalStates3Prongs, std::vector<int> > finalStatesDPlus = 
      {
        {FinalStates3Prongs::KKPi,          std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus}},
        {FinalStates3Prongs::KMinusPiPi,    std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus}},
        {FinalStates3Prongs::KMinusPiPiPi0, std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus, +kPi0}},
        {FinalStates3Prongs::PiPiPi,        std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus}},
        {FinalStates3Prongs::PiPiPiPi0,     std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus, +kPi0}},
      };
    }
    
    // Ds± → K± K∓ π±
    namespace DS {
      enum DecayChannelResoDs {
          PhiPi = 1,
          PhiRho,
          K0starK,
          K0starPi,
          RhoPi,
          RhoK,
          EtaPi,
          f2_1270Pi,
          f2_1370K,
      };
        
      std::unordered_map<DecayChannelResoDs, std::array<int, 2> > resoStatesDs = 
      {
        {DecayChannelResoDs::PhiPi,       std::array<int, 2>{+kPhi,        +kPiPlus}},
        {DecayChannelResoDs::PhiRho,      std::array<int, 2>{+kPhi,        +213}},
        {DecayChannelResoDs::K0starK,     std::array<int, 2>{-kK0Star892,  +kKPlus}},
        {DecayChannelResoDs::K0starPi,    std::array<int, 2>{+kK0Star892,  +kPiPlus}},
        {DecayChannelResoDs::RhoPi,       std::array<int, 2>{113,          +kPiPlus}},
        {DecayChannelResoDs::RhoK,        std::array<int, 2>{113,          +kKPlus}},
        {DecayChannelResoDs::EtaPi,       std::array<int, 2>{221,          +kPiPlus}},
        {DecayChannelResoDs::f2_1270Pi,   std::array<int, 2>{225,          +kPiPlus}},
        {DecayChannelResoDs::f2_1370K,    std::array<int, 2>{10221,        +kKPlus}},
      };

      std::unordered_map<FinalStates3Prongs, std::vector<int> > finalStatesDs = 
      {
        {FinalStates3Prongs::KKPi,       std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus}},
        {FinalStates3Prongs::KKPiPi0,    std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus, +kPi0}},
        {FinalStates3Prongs::KPlusPiPi,  std::vector<int>{+kKPlus,   +kPiPlus,  +kPiMinus}},
        {FinalStates3Prongs::PiPiPi,     std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus}},
        {FinalStates3Prongs::PiPiPiPi0,  std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus, +kPi0}},
      };
    }

    // Dstar → K± K∓ π±
    namespace DStar {
      enum DecayChannelResoDStarD0 {
        RhoPi = 1,
        RhoK,
        K0starPi0,
        K0starPiPlus,
      };

      std::unordered_map<DecayChannelResoDStarD0, std::array<int, 2> > resoStatesDStarD0 = 
      {
        {DecayChannelResoDStarD0::RhoPi,         std::array<int, 2>{213,            +kPiMinus}},
        {DecayChannelResoDStarD0::RhoK,          std::array<int, 2>{213,            +kKMinus}},
        {DecayChannelResoDStarD0::K0starPi0,     std::array<int, 2>{-kK0Star892,    +kPi0}},
        {DecayChannelResoDStarD0::K0starPiPlus,  std::array<int, 2>{-kKPlusStar892, +kPiPlus}},
      };

      std::unordered_map<FinalStates3Prongs, std::vector<int> > finalStatesDStar = 
      {
        {FinalStates3Prongs::KKPi,           std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus}},
        {FinalStates3Prongs::KMinusPiPi,     std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus}},
        {FinalStates3Prongs::KMinusPiPiPi0,  std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus, +kPi0}},
        {FinalStates3Prongs::PiPiPi,         std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus}},
        {FinalStates3Prongs::PiPiPiPi0,      std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus, +kPi0}},
      };
    }

    // Lc → p K∓ π±
    namespace LambdaC {
      enum DecayChannelResoLambdaC {
        K0starP = 1,
        ProtonPhi,
        DeltaK,
        Lambda1520Pi,
      }; 

      std::unordered_map<DecayChannelResoLambdaC, std::array<int, 2> > resoStatesLambdaC = 
      {
        {DecayChannelResoLambdaC::K0starP,        std::array<int, 2>{+kK0Star892, +kProton}},
        {DecayChannelResoLambdaC::ProtonPhi,      std::array<int, 2>{+kProton,    +kPhi}},
        {DecayChannelResoLambdaC::DeltaK,         std::array<int, 2>{+2224,       +kKMinus}},
        {DecayChannelResoLambdaC::Lambda1520Pi,   std::array<int, 2>{+102134,     +kPiPlus}},
      };

      std::unordered_map<FinalStates3Prongs, std::vector<int> > finalStatesLc = 
      {
        {FinalStates3Prongs::ProtonKPi,     std::vector<int>{+kProton,   +kKMinus,  +kPiPlus}},
        {FinalStates3Prongs::ProtonKPiPi0,  std::vector<int>{+kProton,   +kKMinus,  +kPiPlus, +kPi0}},
        {FinalStates3Prongs::ProtonPiPi,    std::vector<int>{+kProton,   +kPiMinus, +kPiPlus}},
        {FinalStates3Prongs::ProtonKK,      std::vector<int>{+kProton,   +kKMinus,  +kKPlus}}
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
        {DecayChannelResoXiC::K0starP,    std::array<int, 2>{-kK0Star892, +kProton}},
        {DecayChannelResoXiC::ProtonPhi,  std::array<int, 2>{+kProton,    +kPhi}},
        {DecayChannelResoXiC::SigmaPiPi,  std::array<int, 2>{+kSigmaPlus, +kPiPlus}},
      };

      std::unordered_map<FinalStates3Prongs, std::vector<int> > finalStatesXic = 
      {
        {FinalStates3Prongs::ProtonKPi,  std::vector<int>{+kProton,   +kKMinus,  +kPiPlus}},
        {FinalStates3Prongs::ProtonKK,   std::vector<int>{+kProton,   +kKMinus,  +kKPlus}},
        {FinalStates3Prongs::SigmaPiPi,  std::vector<int>{+kSigmaPlus,   +kPiMinus,  +kPiPlus}},
      };
    }

    std::unordered_map<FinalStates3Prongs, std::vector<int> > getParticleFinalStates3Prongs(int pdgMother) {
      switch (pdgMother) {
        case Pdg::kDPlus:
          return DPlus::finalStatesDPlus;
        case Pdg::kDS:
          return DS::finalStatesDs;
        case Pdg::kDStar:
          return DStar::finalStatesDStar;
        case Pdg::kLambdaCPlus:
          return LambdaC::finalStatesLc;
        case Pdg::kXiCPlus:
          return XiC::finalStatesXic;
        default:
          LOG(error) << "Unknown PDG code for 3-prong final states: " << pdgMother;
          return {};
      }
    }

    bool checkResonantDecay(std::vector<int> arrDaughIndex, std::array<int, 2> arrPDGResonant) {
      // LOG(info) << "Entered checkResonantDecay with daughters: " << arrDaughIndex[0] << ", " << arrDaughIndex[1] << " and resonant PDG codes: " << arrPDGResonant[0] << ", " << arrPDGResonant[1];
      // LOG(info) << "arrDaughIndex.size(): " << arrDaughIndex.size() << ", arrPDGResonant.size(): " << arrPDGResonant.size();
      std::array<int, 2> arrPDGResonantAbs = {std::abs(arrPDGResonant[0]), std::abs(arrPDGResonant[1])};
      LOG(info) << "Testing: " << arrDaughIndex[0] << ", " << arrDaughIndex[1] << " matching PDG codes: " << arrPDGResonant[0] << ", " << arrPDGResonant[1];
      for (int i = 0; i < 2; i++) {
        LOG(info) << "Checking daughter index: " << arrDaughIndex[i];
        bool findDaug = false;
        for (int j = 0; j < 2; j++) {
          LOG(info) << "Checking daughter PDG: " << arrDaughIndex[i] << " against resonant PDG: " << arrPDGResonantAbs[j];
          if (std::abs(arrDaughIndex[i]) == arrPDGResonantAbs[j]) {
            arrPDGResonantAbs[j] = -1; // Mark as found
            LOG(info) << "Matched!";
            findDaug = true;
            break; // If a daughter matches, break the inner loop
          }
        }
        if (!findDaug) {
          LOG(info) << "Returning false";
          return false; // If any daughter does not match, return false
        }
      }
      LOG(info) << "Resonant decay found with daughters: " << arrDaughIndex[0] << ", " << arrDaughIndex[1] << " matching PDG codes: " << arrPDGResonant[0] << ", " << arrPDGResonant[1];
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
            std::cout << "Checking D0 resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << " vs " << arrDaughIndex[0] << " " << arrDaughIndex[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "D0 resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kDPlus:
          for (const auto& [flag, pdgCodes] : DPlus::resoStatesDPlus) {
            // std::cout << "Checking DPlus resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "D+ resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kDS:
          for (const auto& [flag, pdgCodes] : DS::resoStatesDs) {
            // std::cout << "Checking DS resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "Ds resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kDStar:
          for (const auto& [flag, pdgCodes] : DStar::resoStatesDStarD0) {
            std::cout << "Checking DStar resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "Dstar resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kLambdaCPlus:
          for (const auto& [flag, pdgCodes] : LambdaC::resoStatesLambdaC) {
            // std::cout << "Checking LambdaC resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "Lc resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kXiCPlus:
          for (const auto& [flag, pdgCodes] : XiC::resoStatesXiC) {
            // std::cout << "Checking XiC resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "Xic resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
      }
      LOG(info) << "Leaving function with channel: " << static_cast<int>(*channel);
    }

} // namespace o2::hf_corrbkg

#endif // PWGHF_CORE_CORRELATEDBKGS_H_
