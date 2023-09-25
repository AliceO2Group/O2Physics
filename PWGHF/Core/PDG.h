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
constexpr double MassB0 = 5.27953;
constexpr double MassB0Bar = 5.27953;
constexpr double MassBPlus = 5.27915;
constexpr double MassBS = 5.3663;
constexpr double MassBSBar = 5.3663;
constexpr double MassChiC1 = 3.51066;
constexpr double MassD0 = 1.86484;
constexpr double MassD0Bar = 1.86484;
constexpr double MassDMinus = 1.86962;
constexpr double MassDPlus = 1.86962;
constexpr double MassDS = 1.9685;
constexpr double MassDSBar = 1.9685;
constexpr double MassDStar = 2.01027;
constexpr double MassJPsi = 3.09692;
constexpr double MassLambdaB0 = 5.6202;
constexpr double MassLambdaCPlus = 2.28646;
constexpr double MassOmegaC0 = 2.6975;
constexpr double MassPhi = 1.01946;
constexpr double MassSigmaC0 = 2.45376;
constexpr double MassSigmaCPlusPlus = 2.45402;
constexpr double MassX3872 = 3.87165;
constexpr double MassXiCCPlusPlus = 3.59798;
constexpr double MassXiCPlus = 2.4679;
constexpr double MassXiCZero = 2.471;

// constexpr double massPi = o2::constants::physics::MassPionCharged;
// constexpr double massKa = o2::constants::physics::MassKaonCharged;
// constexpr double massProton = o2::constants::physics::MassProton;
// constexpr double massGamma = o2::constants::physics::MassPhoton;
// constexpr double massK0S = o2::constants::physics::MassKaonNeutral;
// constexpr double massLambda = o2::constants::physics::MassLambda;
// constexpr double massXi = o2::constants::physics::MassXiMinus;
// constexpr double massPhi = 1.019455;
// constexpr double massD0 = 1.86484;
// constexpr double massDPlus = 1.86962;
// constexpr double massDS = 1.9685;
// constexpr double massLambdaCPlus = 2.28646;
// constexpr double massXiCPlus = 2.4679;
// constexpr double massDStar = 2.01027;
// constexpr double massBPlus = 5.27915;
// constexpr double massB0 = 5.27953;
// constexpr double massBs = 5.3663;
// constexpr double massLb = 5.6202;
// constexpr double massXib = 5.7924;
} // namespace o2::analysis::pdg

#endif // PWGHF_CORE_PDG_H_
