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
/// \file   mcParticle.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \since  31/05/2024
/// \brief  Utilities to handle the MC information
///

#ifndef PWGLF_UTILS_MCPARTICLE_H_
#define PWGLF_UTILS_MCPARTICLE_H_

#include "ReconstructionDataFormats/PID.h"

#include <string>
#include <vector>

namespace o2
{
namespace pwglf
{

/// @brief  Extended the indices of the PID class with additional particles
class PIDExtended
{
 public:
  typedef int16_t ID;

  static constexpr ID Electron = 0;
  static constexpr ID Muon = 1;
  static constexpr ID Pion = 2;
  static constexpr ID Kaon = 3;
  static constexpr ID Proton = 4;
  static constexpr ID Deuteron = 5;
  static constexpr ID Triton = 6;
  static constexpr ID Helium3 = 7;
  static constexpr ID Alpha = 8;
  static constexpr ID PI0 = 9;
  static constexpr ID Photon = 10;
  static constexpr ID K0 = 11;
  static constexpr ID Lambda = 12;
  static constexpr ID HyperTriton = 13;
  static constexpr ID Hyperhydrog4 = 14;
  static constexpr ID XiMinus = 15;
  static constexpr ID OmegaMinus = 16;
  static constexpr ID HyperHelium4 = 17;
  static constexpr ID HyperHelium5 = 18;

  static_assert(Electron == o2::track::PID::Electron, "PID::Electron mismatch");
  static_assert(Muon == o2::track::PID::Muon, "PID::Muon mismatch");
  static_assert(Pion == o2::track::PID::Pion, "PID::Pion mismatch");
  static_assert(Kaon == o2::track::PID::Kaon, "PID::Kaon mismatch");
  static_assert(Proton == o2::track::PID::Proton, "PID::Proton mismatch");
  static_assert(Deuteron == o2::track::PID::Deuteron, "PID::Deuteron mismatch");
  static_assert(Triton == o2::track::PID::Triton, "PID::Triton mismatch");
  static_assert(Helium3 == o2::track::PID::Helium3, "PID::Helium3 mismatch");
  static_assert(Alpha == o2::track::PID::Alpha, "PID::Alpha mismatch");
  static_assert(PI0 == o2::track::PID::PI0, "PID::PI0 mismatch");
  static_assert(Photon == o2::track::PID::Photon, "PID::Photon mismatch");
  static_assert(K0 == o2::track::PID::K0, "PID::K0 mismatch");
  static_assert(Lambda == o2::track::PID::Lambda, "PID::Lambda mismatch");
  static_assert(HyperTriton == o2::track::PID::HyperTriton, "PID::HyperTriton mismatch");
  static_assert(Hyperhydrog4 == o2::track::PID::Hyperhydrog4, "PID::Hyperhydrog4 mismatch");
  static_assert(XiMinus == o2::track::PID::XiMinus, "PID::XiMinus mismatch");
  static_assert(OmegaMinus == o2::track::PID::OmegaMinus, "PID::OmegaMinus mismatch");
  // static_assert(HyperHelium4 == o2::track::PID::HyperHelium4, "PID::HyperHelium4 mismatch");
  // static_assert(HyperHelium5 == o2::track::PID::HyperHelium5, "PID::HyperHelium5 mismatch");

  static constexpr ID PIDCountsUntilAl = 9; // Number of indices defined in PID.h equivalent to o2::track::PID::NIDs
  // static_assert(PIDCountsUntilAl == o2::track::PID::NIDs, "PID::NIDs mismatch");

  static constexpr ID PIDCounts = 19; // Number of indices defined in PID.h
  // static_assert(PIDCounts == o2::track::PID::NIDsTot, "PID::NIDsTot mismatch");
  // Define an array of IDs
  static constexpr std::array<ID, PIDCountsUntilAl> mIDsUntilAl = {Electron, Muon, Pion, Kaon, Proton, Deuteron, Triton, Helium3, Alpha};
  static constexpr std::array<ID, PIDCounts> mIDs = {Electron, Muon, Pion, Kaon, Proton, Deuteron, Triton, Helium3, Alpha, PI0, Photon, K0, Lambda, HyperTriton, Hyperhydrog4, XiMinus, OmegaMinus, HyperHelium4, HyperHelium5};

  // Define the antiparticles
  static constexpr ID Positron = PIDCounts;
  static constexpr ID MuonPlus = PIDCounts + 1;
  static constexpr ID PionMinus = PIDCounts + 2;
  static constexpr ID KaonMinus = PIDCounts + 3;
  static constexpr ID AntiProton = PIDCounts + 4;
  static constexpr ID AntiDeuteron = PIDCounts + 5;
  static constexpr ID AntiTriton = PIDCounts + 6;
  static constexpr ID AntiHelium3 = PIDCounts + 7;
  static constexpr ID AntiAlpha = PIDCounts + 8;
  static constexpr ID AntiLambda = PIDCounts + 9;
  static constexpr ID AntiHyperTriton = PIDCounts + 10;
  static constexpr ID AntiHyperhydrog4 = PIDCounts + 11;
  static constexpr ID XiPlus = PIDCounts + 12;
  static constexpr ID OmegaPlus = PIDCounts + 13;
  static constexpr ID AntiHyperHelium4 = PIDCounts + 14;
  static constexpr ID AntiHyperHelium5 = PIDCounts + 15;

  static constexpr ID Neutron = PIDCounts + 16;
  static constexpr ID AntiNeutron = PIDCounts + 17;
  static constexpr ID Phi = PIDCounts + 18;
  static constexpr ID BZero = PIDCounts + 19;
  static constexpr ID BPlus = PIDCounts + 20;
  static constexpr ID BS = PIDCounts + 21;
  static constexpr ID D0 = PIDCounts + 22;
  static constexpr ID DPlus = PIDCounts + 23;
  static constexpr ID DS = PIDCounts + 24;
  static constexpr ID DStar = PIDCounts + 25;
  static constexpr ID ChiC1 = PIDCounts + 26;
  static constexpr ID JPsi = PIDCounts + 27;
  static constexpr ID LambdaB0 = PIDCounts + 28;
  static constexpr ID LambdaCPlus = PIDCounts + 29;
  static constexpr ID OmegaC0 = PIDCounts + 30;
  static constexpr ID SigmaC0 = PIDCounts + 31;
  static constexpr ID SigmaCPlusPlus = PIDCounts + 32;
  static constexpr ID X3872 = PIDCounts + 33;
  static constexpr ID Xi0 = PIDCounts + 34;
  static constexpr ID XiB0 = PIDCounts + 35;
  static constexpr ID XiCCPlusPlus = PIDCounts + 36;
  static constexpr ID XiCPlus = PIDCounts + 37;
  static constexpr ID XiC0 = PIDCounts + 38;
  static constexpr ID NIDsTot = PIDCounts + 39;

  static constexpr const char* sNames[NIDsTot + 1] = {
    o2::track::pid_constants::sNames[Electron],     // Electron
    o2::track::pid_constants::sNames[Muon],         // Muon
    o2::track::pid_constants::sNames[Pion],         // Pion
    o2::track::pid_constants::sNames[Kaon],         // Kaon
    o2::track::pid_constants::sNames[Proton],       // Proton
    o2::track::pid_constants::sNames[Deuteron],     // Deuteron
    o2::track::pid_constants::sNames[Triton],       // Triton
    o2::track::pid_constants::sNames[Helium3],      // Helium3
    o2::track::pid_constants::sNames[Alpha],        // Alpha
    o2::track::pid_constants::sNames[PI0],          // PI0
    o2::track::pid_constants::sNames[Photon],       // Photon
    o2::track::pid_constants::sNames[K0],           // K0
    o2::track::pid_constants::sNames[Lambda],       // Lambda
    o2::track::pid_constants::sNames[HyperTriton],  // HyperTriton
    o2::track::pid_constants::sNames[Hyperhydrog4], // Hyperhydrog4
    o2::track::pid_constants::sNames[XiMinus],      // XiMinus
    o2::track::pid_constants::sNames[OmegaMinus],   // OmegaMinus
    "HyperHelium4",                                 // HyperHelium4
    "HyperHelium5",                                 // HyperHelium5
    "Positron",                                     // Positron
    "MuonPlus",                                     // MuonPlus
    "PionMinus",                                    // PionMinus
    "KaonMinus",                                    // KaonMinus
    "AntiProton",                                   // AntiProton
    "AntiDeuteron",                                 // AntiDeuteron
    "AntiTriton",                                   // AntiTriton
    "AntiHelium3",                                  // AntiHelium3
    "AntiAlpha",                                    // AntiAlpha
    "AntiLambda",                                   // AntiLambda
    "AntiHyperTriton",                              // AntiHyperTriton
    "AntiHyperhydrog4",                             // AntiHyperhydrog4
    "XiPlus",                                       // XiPlus
    "OmegaPlus",                                    // OmegaPlus
    "AntiHyperHelium4",                             // AntiHyperHelium4
    "AntiHyperHelium5",                             // AntiHyperHelium5
    "Neutron",                                      // Neutron
    "AntiNeutron",                                  // AntiNeutron
    "Phi",                                          // Phi
    "BZero",                                        // BZero
    "BPlus",                                        // BPlus
    "BS",                                           // BS
    "D0",                                           // D0
    "DPlus",                                        // DPlus
    "DS",                                           // DS
    "DStar",                                        // DStar
    "ChiC1",                                        // ChiC1
    "JPsi",                                         // JPsi
    "LambdaB0",                                     // LambdaB0
    "LambdaCPlus",                                  // LambdaCPlus
    "OmegaC0",                                      // OmegaC0
    "SigmaC0",                                      // SigmaC0
    "SigmaCPlusPlus",                               // SigmaCPlusPlus
    "X3872",                                        // X3872
    "Xi0",                                          // Xi0
    "XiB0",                                         // XiB0
    "XiCCPlusPlus",                                 // XiCCPlusPlus
    "XiCPlus",                                      // XiCPlus
    "XiC0",                                         // XiC0
    nullptr};

  static std::vector<std::string> arrayNames()
  {
    std::vector<std::string> names;
    for (int i = 0; i < NIDsTot; i++) {
      names.push_back(sNames[i]);
    }
    return names;
  }

  static const char* getName(ID id) { return sNames[id]; }

  /// \brief Convert PDG code to PID index
  template <typename TrackType>
  static ID pdgToId(const TrackType& particle)
  {
    switch (particle.pdgCode()) {
      case 11:
        return Electron;
      case -11:
        return Positron;
      case 13:
        return Muon;
      case -13:
        return MuonPlus;
      case 211:
        return Pion;
      case -211:
        return PionMinus;
      case 321:
        return Kaon;
      case -321:
        return KaonMinus;
      case 2212:
        return Proton;
      case -2212:
        return AntiProton;
      case 2112:
        return Neutron;
      case -2112:
        return AntiNeutron;
      case o2::constants::physics::Pdg::kDeuteron:
        return Deuteron;
      case -o2::constants::physics::Pdg::kDeuteron:
        return AntiDeuteron;
      case o2::constants::physics::Pdg::kTriton:
        return Triton;
      case -o2::constants::physics::Pdg::kTriton:
        return AntiTriton;
      case o2::constants::physics::Pdg::kHelium3:
        return Helium3;
      case -o2::constants::physics::Pdg::kHelium3:
        return AntiHelium3;
      case o2::constants::physics::Pdg::kAlpha:
        return Alpha;
      case -o2::constants::physics::Pdg::kAlpha:
        return AntiAlpha;
      case o2::constants::physics::Pdg::kHyperHelium4:
        return HyperHelium4;
      case -o2::constants::physics::Pdg::kHyperHelium4:
        return AntiHyperHelium4;
      case o2::constants::physics::Pdg::kHyperHelium5:
        return HyperHelium5;
      case -o2::constants::physics::Pdg::kHyperHelium5:
        return AntiHyperHelium5;
      case 111:
        return PI0;
      case 22:
        return Photon;
      case 130:
        return K0;
      case 3122:
        return Lambda;
      case -3122:
        return AntiLambda;
      case o2::constants::physics::Pdg::kHyperTriton:
        return HyperTriton;
      case -o2::constants::physics::Pdg::kHyperTriton:
        return AntiHyperTriton;
      case o2::constants::physics::Pdg::kHyperHydrogen4:
        return Hyperhydrog4;
      case -o2::constants::physics::Pdg::kHyperHydrogen4:
        return AntiHyperhydrog4;
      case 3312:
        return XiMinus;
      case -3312:
        return XiPlus;
      case 3334:
        return OmegaMinus;
      case -3334:
        return OmegaPlus;
      case o2::constants::physics::Pdg::kPhi:
        return Phi;
      case o2::constants::physics::Pdg::kB0:
        return BZero;
      case o2::constants::physics::Pdg::kBPlus:
        return BPlus;
      case o2::constants::physics::Pdg::kBS:
        return BS;
      case o2::constants::physics::Pdg::kD0:
        return D0;
      case o2::constants::physics::Pdg::kDPlus:
        return DPlus;
      case o2::constants::physics::Pdg::kDS:
        return DS;
      case o2::constants::physics::Pdg::kDStar:
        return DStar;
      case o2::constants::physics::Pdg::kChiC1:
        return ChiC1;
      case o2::constants::physics::Pdg::kJPsi:
        return JPsi;
      case o2::constants::physics::Pdg::kLambdaB0:
        return LambdaB0;
      case o2::constants::physics::Pdg::kLambdaCPlus:
        return LambdaCPlus;
      case o2::constants::physics::Pdg::kOmegaC0:
        return OmegaC0;
      case o2::constants::physics::Pdg::kSigmaC0:
        return SigmaC0;
      case o2::constants::physics::Pdg::kSigmaCPlusPlus:
        return SigmaCPlusPlus;
      case o2::constants::physics::Pdg::kX3872:
        return X3872;
      case o2::constants::physics::Pdg::kXi0:
        return Xi0;
      case o2::constants::physics::Pdg::kXiB0:
        return XiB0;
      case o2::constants::physics::Pdg::kXiCCPlusPlus:
        return XiCCPlusPlus;
      case o2::constants::physics::Pdg::kXiCPlus:
        return XiCPlus;
      case o2::constants::physics::Pdg::kXiC0:
        return XiC0;
      default:
        LOG(debug) << "Cannot identify particle with PDG code " << particle.pdgCode();
        break;
    }
    return -1;
  }
};

} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_MCPARTICLE_H_
