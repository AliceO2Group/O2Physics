// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUtils.h
/// \brief Collision selection
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_FEMTOUTILS_H_
#define PWGCF_FEMTO_CORE_FEMTOUTILS_H_

#include "RecoDecay.h"

#include "Common/Core/TableHelper.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/InitContext.h"

#include "TPDGCode.h"

#include "fairlogger/Logger.h"

#include <cmath>
#include <cstdint>
#include <experimental/type_traits>
#include <optional>
#include <unordered_map>
#include <utility>

namespace o2::analysis::femto
{
namespace utils
{

template <typename T1, typename T2>
inline std::optional<T2> getIndex(const T1& index, const std::unordered_map<T1, T2>& map)
{
  auto it = map.find(index);
  if (it != map.end()) {
    return it->second;
  }
  return std::nullopt;
}

template <typename T>
float itsSignal(T const& track)
{
  uint32_t clsizeflag = track.itsClusterSizes();
  auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
  auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
  auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
  auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
  auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
  auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
  auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;
  int numLayers = 7;
  int sumClusterSizes = clSizeLayer0 + clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6;
  float cosLamnda = 1. / std::cosh(track.eta());
  return (static_cast<float>(sumClusterSizes) / numLayers) * cosLamnda;
};

template <typename T>
float sphericity(T const& tracks)
{

  int minNumberTracks = 2;
  float maxSphericity = 2.f;

  if (tracks.size() <= minNumberTracks) {
    return maxSphericity;
  }

  // Initialize the transverse momentum tensor components
  float sxx = 0.f;
  float syy = 0.f;
  float sxy = 0.f;
  float sumPt = 0.f;

  // Loop over the tracks to compute the tensor components
  for (const auto& track : tracks) {
    sxx += (track.px() * track.px()) / track.pt();
    syy += (track.py() * track.py()) / track.pt();
    sxy += (track.px() * track.py()) / track.pt();
    sumPt += track.pt();
  }
  sxx /= sumPt;
  syy /= sumPt;
  sxy /= sumPt;

  // Compute the eigenvalues (real values)
  float lambda1 = ((sxx + syy) + std::sqrt((sxx + syy) * (sxx + syy) - 4 * (sxx * syy - sxy * sxy))) / 2;
  float lambda2 = ((sxx + syy) - std::sqrt((sxx + syy) * (sxx + syy) - 4 * (sxx * syy - sxy * sxy))) / 2;

  if (lambda1 <= 0.f || lambda2 <= 0.f) {
    return maxSphericity;
  }

  // Compute sphericity
  return 2.f * lambda2 / (lambda1 + lambda2);
}

inline float getMass(int pdgCode)
{
  // use this function instead of TDatabasePDG to return masses defined in the PhysicsConstants.h header
  // this approach saves a lot of memory and important partilces like deuteron are missing in TDatabasePDG anyway
  float mass = 0.f;
  // add new particles if necessary here
  switch (std::abs(pdgCode)) {
    case kPiPlus:
      mass = o2::constants::physics::MassPiPlus;
      break;
    case kKPlus:
      mass = o2::constants::physics::MassKPlus;
      break;
    case kProton:
      mass = o2::constants::physics::MassProton;
      break;
    case kLambda0:
      mass = o2::constants::physics::MassLambda;
      break;
    case o2::constants::physics::Pdg::kPhi:
      mass = o2::constants::physics::MassPhi;
      break;
    case kRho770_0:
      mass = 775.26; // not defined in O2?
      break;
    case kRho770Plus:
      mass = 775.11; // not defined in O2?
      break;
    case o2::constants::physics::Pdg::kK0Star892:
      mass = o2::constants::physics::MassK0Star892;
      break;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      mass = o2::constants::physics::MassLambdaCPlus;
      break;
    case o2::constants::physics::Pdg::kDeuteron:
      mass = o2::constants::physics::MassDeuteron;
      break;
    case o2::constants::physics::Pdg::kTriton:
      mass = o2::constants::physics::MassTriton;
      break;
    case o2::constants::physics::Pdg::kHelium3:
      mass = o2::constants::physics::MassHelium3;
      break;
    case kSigmaMinus:
      mass = o2::constants::physics::MassSigmaMinus;
      break;
    case kSigmaPlus:
      mass = o2::constants::physics::MassSigmaPlus;
      break;
    case kXiMinus:
      mass = o2::constants::physics::MassXiMinus;
      break;
    case kOmegaMinus:
      mass = o2::constants::physics::MassOmegaMinus;
      break;
    default:
      LOG(fatal) << "PDG code is not suppored";
  }
  return mass;
}

template <typename T>
float qn(T const& col)
{
  float qn = std::sqrt(col.qvecFT0CReVec()[0] * col.qvecFT0CReVec()[0] + col.qvecFT0CImVec()[0] * col.qvecFT0CImVec()[0]) * std::sqrt(col.sumAmplFT0C());
  return qn;
}

inline std::optional<float> dphistar(float magfield, float radius, float signedPt, float phi)
{
  float arg = 0.3f * (0.1f * magfield) * (0.01 * radius) / (2.f * signedPt);
  if (std::fabs(arg) <= 1.f) {
    return RecoDecay::constrainAngle(phi - std::asin(arg));
  }
  return std::nullopt;
}

inline bool enableTable(const char* tableName, int userSetting, o2::framework::InitContext& initContext)
{
  if (userSetting == 1) {
    LOG(info) << "Enabled femto table (forced on): " << tableName;
    return true;
  }
  if (userSetting == 0) {
    LOG(info) << "Disabled femto table (forced off): " << tableName;
    return false;
  }
  bool required = o2::common::core::isTableRequiredInWorkflow(initContext, tableName);
  if (required) {
    LOG(info) << "Enabled femto table (auto): " << tableName;
  }
  return required;
}

template <typename T>
using HasMass = decltype(std::declval<T&>().mass());

template <typename T>
using HasSign = decltype(std::declval<T&>().sign());

}; // namespace utils
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_FEMTOUTILS_H_
