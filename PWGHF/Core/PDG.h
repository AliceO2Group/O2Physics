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

/// \file PDG.h
/// \brief PDG features for HF

#ifndef PWGHF_CORE_PDG_H_
#define PWGHF_CORE_PDG_H_

namespace o2::analysis::pdg
{
/// \brief Declarations of named PDG codes of HF particles missing in ROOT PDG_t.
/// \note Follow kCamelCase naming convention
/// \link https://root.cern/doc/master/TPDGCode_8h.html
enum Code {
  kB0 = 511,
  kB0Bar = -511,
  kBPlus = 521,
  kBS = 531,
  kBSBar = -531,
  kD0 = 421,
  kD0Bar = -421,
  kDMinus = -411,
  kDPlus = 411,
  kDS = 431,
  kDSBar = -431,
  kDStar = 413,
  kChiC1 = 20443,
  kJPsi = 443,
  kLambdaB0 = 5122,
  kLambdaCPlus = 4122,
  kOmegaC0 = 4332,
  kPhi = 333,
  kSigmaC0 = 4112,
  kSigmaCPlusPlus = 4222,
  kX3872 = 9920443,
  kXiCCPlusPlus = 4422,
  kXiCPlus = 4232,
  kXiCZero = 4132
};

/// \brief Declarations of particle massed from o2::O2DatabasePDG
/// \note Generated with root -b -l -q PWGHF/Core/listMasses.C
constexpr float MassB0 = 5.27953;
constexpr float MassB0Bar = 5.27953;
constexpr float MassBPlus = 5.27915;
constexpr float MassBS = 5.3663;
constexpr float MassBSBar = 5.3663;
constexpr float MassChiC1 = 3.51066;
constexpr float MassD0 = 1.86484;
constexpr float MassD0Bar = 1.86484;
constexpr float MassDMinus = 1.86962;
constexpr float MassDPlus = 1.86962;
constexpr float MassDS = 1.9685;
constexpr float MassDSBar = 1.9685;
constexpr float MassDStar = 2.01027;
constexpr float MassJPsi = 3.09692;
constexpr float MassLambdaB0 = 5.6202;
constexpr float MassLambdaCPlus = 2.28646;
constexpr float MassOmegaC0 = 2.6975;
constexpr float MassPhi = 1.01946;
constexpr float MassSigmaC0 = 2.45376;
constexpr float MassSigmaCPlusPlus = 2.45402;
constexpr float MassX3872 = 3.87165;
constexpr float MassXiCCPlusPlus = 3.59798;
constexpr float MassXiCPlus = 2.4679;
constexpr float MassXiCZero = 2.471;

// constexpr float massPi = o2::constants::physics::MassPionCharged;
// constexpr float massKa = o2::constants::physics::MassKaonCharged;
// constexpr float massProton = o2::constants::physics::MassProton;
// constexpr float massGamma = o2::constants::physics::MassPhoton;
// constexpr float massK0S = o2::constants::physics::MassKaonNeutral;
// constexpr float massLambda = o2::constants::physics::MassLambda;
// constexpr float massXi = o2::constants::physics::MassXiMinus;
// constexpr float massPhi = 1.019455;
// constexpr float massD0 = 1.86484;
// constexpr float massDPlus = 1.86962;
// constexpr float massDS = 1.9685;
// constexpr float massLambdaCPlus = 2.28646;
// constexpr float massXiCPlus = 2.4679;
// constexpr float massDStar = 2.01027;
// constexpr float massBPlus = 5.27915;
// constexpr float massB0 = 5.27953;
// constexpr float massBs = 5.3663;
// constexpr float massLb = 5.6202;
// constexpr float massXib = 5.7924;
} // namespace o2::analysis::pdg

#endif // PWGHF_CORE_PDG_H_
