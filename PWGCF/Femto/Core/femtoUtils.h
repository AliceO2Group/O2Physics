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

#include "Common/Core/TableHelper.h"

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
  double cosLamnda = 1. / std::cosh(track.eta());
  double signal = (static_cast<double>(sumClusterSizes) / numLayers) * cosLamnda;
  return static_cast<float>(signal);
};

inline double getMass(int pdgCode)
{
  // use this function instead of TDatabasePDG to return masses defined in the PhysicsConstants.h header
  // this approach saves a lot of memory and important partilces like deuteron are missing in TDatabasePDG anyway
  double mass = 0.f;
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

/// Recalculate pT for Kinks (Sigmas) using kinematic constraints
inline float calcPtnew(float pxMother, float pyMother, float pzMother, float pxDaughter, float pyDaughter, float pzDaughter)
{
  // Particle masses in GeV/c^2
  const float massPion = 0.13957f;
  const float massNeutron = 0.93957f;
  const float massSigmaMinus = 1.19745f;

  // Calculate mother momentum and direction versor
  float pMother = std::sqrt(pxMother * pxMother + pyMother * pyMother + pzMother * pzMother);
  if (pMother < 1e-6f)
    return -999.f;

  float versorX = pxMother / pMother;
  float versorY = pyMother / pMother;
  float versorZ = pzMother / pMother;

  // Calculate daughter energy
  float ePi = std::sqrt(massPion * massPion + pxDaughter * pxDaughter + pyDaughter * pyDaughter + pzDaughter * pzDaughter);

  // Scalar product of versor with daughter momentum
  float a = versorX * pxDaughter + versorY * pyDaughter + versorZ * pzDaughter;

  // Solve quadratic equation for momentum magnitude
  float K = massSigmaMinus * massSigmaMinus + massPion * massPion - massNeutron * massNeutron;
  float A = 4.f * (ePi * ePi - a * a);
  float B = -4.f * a * K;
  float C = 4.f * ePi * ePi * massSigmaMinus * massSigmaMinus - K * K;

  if (std::abs(A) < 1e-6f)
    return -999.f;

  float D = B * B - 4.f * A * C;
  if (D < 0.f)
    return -999.f;

  float sqrtD = std::sqrt(D);
  float P1 = (-B + sqrtD) / (2.f * A);
  float P2 = (-B - sqrtD) / (2.f * A);

  // Pick physical solution: prefer P2 if positive, otherwise P1
  if (P2 < 0.f && P1 < 0.f)
    return -999.f;
  if (P2 < 0.f)
    return P1;

  // Choose solution closest to original momentum
  float p1Diff = std::abs(P1 - pMother);
  float p2Diff = std::abs(P2 - pMother);
  float P = (p1Diff < p2Diff) ? P1 : P2;

  // Calculate pT from recalibrated momentum
  float pxS = versorX * P;
  float pyS = versorY * P;
  return std::sqrt(pxS * pxS + pyS * pyS);
}

/// Helper function to calculate recalculated pT for kink particles (Sigma/SigmaPlus)
template <typename TParticle, typename TTrackTable>
float getRecalculatedPtForKink(const TParticle& particle, const TTrackTable& trackTable)
{
  // Check if particle has chaDau index
  if constexpr (requires { particle.has_chaDau(); }) {
    if (particle.has_chaDau()) {
      try {
        auto chaDaughter = trackTable.rawIteratorAt(particle.chaDauId() - trackTable.offset());

        // Extract momentum components directly from dynamic columns
        float pxDaug = chaDaughter.px();
        float pyDaug = chaDaughter.py();
        float pzDaug = chaDaughter.pz();

        // Get momentum components from dynamic columns
        float pxMoth = particle.px();
        float pyMoth = particle.py();
        float pzMoth = particle.pz();

        // Recalculate pT using kinematic constraints
        float ptRecalc = calcPtnew(pxMoth, pyMoth, pzMoth, pxDaug, pyDaug, pzDaug);
        if (ptRecalc > 0) {
          return ptRecalc;
        }
      } catch (const std::exception& e) {
        return -1.0f;
      }
    }
  }
  return -1.0f;
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

template <typename T>
inline int signum(T x)
{
  return (T(0) < x) - (x < T(0));
}

}; // namespace utils
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_FEMTOUTILS_H_
