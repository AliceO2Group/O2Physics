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
/// \file   Decayer.h
/// \author Jesper Karlsson Gumprecht
/// \since  15/12/2025
/// \brief  Basic class to handle short-lived particle decays in the fast simulation
///

#ifndef ALICE3_CORE_DECAYER_H_
#define ALICE3_CORE_DECAYER_H_

#include "ALICE3/Core/OTFParticle.h"
#include "ALICE3/Core/TrackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/Track.h>

#include <TDecayChannel.h> // IWYU pragma: keep
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace o2::upgrade
{

class Decayer
{
 public:
  // Default constructor
  Decayer() = default;

  template <typename TDatabase>
  std::vector<o2::upgrade::OTFParticle> decayParticle(const OTFParticle& particle, const TDatabase& pdgDB)
  {
    auto particleInfo = pdgDB->GetParticle(particle.pdgCode());
    if (!particleInfo) {
      return {};
    }

    const int charge = particleInfo->Charge() / 3;
    const double mass = particleInfo->Mass();
    std::array<double, 3> decayVtx = generateDecayVertex<double>(particle, pdgDB);
    mVx = decayVtx[0];
    mVy = decayVtx[1];
    mVz = decayVtx[2];
    double px{}, py{}, e{};

    if (!charge) {
      px = particle.px();
      py = particle.py();
    } else {
      px = particle.px() * std::cos(mTheta) - particle.py() * std::sin(mTheta);
      py = particle.py() * std::cos(mTheta) + particle.px() * std::sin(mTheta);
    }

    double brTotal = 0.;
    e = std::sqrt(mass * mass + px * px + py * py + particle.pz() * particle.pz());
    for (int ch = 0; ch < particleInfo->NDecayChannels(); ++ch) {
      brTotal += particleInfo->DecayChannel(ch)->BranchingRatio();
    }

    double brSum = 0.;
    std::vector<double> dauMasses;
    std::vector<int> pdgCodesDaughters;
    const double randomChannel = mRand3.Uniform(0., brTotal);
    for (int ch = 0; ch < particleInfo->NDecayChannels(); ++ch) {
      brSum += particleInfo->DecayChannel(ch)->BranchingRatio();
      if (randomChannel < brSum) {
        for (int dau = 0; dau < particleInfo->DecayChannel(ch)->NDaughters(); ++dau) {
          const int pdgDau = particleInfo->DecayChannel(ch)->DaughterPdgCode(dau);
          pdgCodesDaughters.push_back(pdgDau);
          const auto& dauInfo = pdgDB->GetParticle(pdgDau);
          dauMasses.push_back(dauInfo->Mass());
        }
        break;
      }
    }

    if (dauMasses.empty()) {
      return {};
    }

    TLorentzVector tlv(px, py, particle.pz(), e);
    TGenPhaseSpace decay;
    decay.SetDecay(tlv, dauMasses.size(), dauMasses.data());
    decay.Generate();

    std::vector<o2::upgrade::OTFParticle> decayProducts;
    for (size_t i = 0; i < dauMasses.size(); ++i) {
      o2::upgrade::OTFParticle particle;
      TLorentzVector dau = *decay.GetDecay(i);
      particle.setPDG(pdgCodesDaughters[i]);
      particle.setVxVyVz(mVx, mVy, mVz);
      particle.setPxPyPzE(dau.Px(), dau.Py(), dau.Pz(), dau.E());
      particle.setBitOn(o2::upgrade::DecayerBits::ProducedByDecayer);
      decayProducts.push_back(particle);
    }

    return decayProducts;
  }

  template <typename T = float, typename TDatabase, typename TParticle>
  std::array<T, 3> generateDecayVertex(const TParticle& particle, const TDatabase& pdgDB)
  {
    std::array<T, 3> decayVertex{};
    auto particleInfo = pdgDB->GetParticle(particle.pdgCode());
    if (!particleInfo) {
      return {};
    }

    const int charge = particleInfo->Charge() / 3;
    const double mass = particleInfo->Mass();
    const double u = mRand3.Uniform(0.001, 0.999);
    const double ctau = o2::constants::physics::LightSpeedCm2S * particleInfo->Lifetime(); // cm
    const double betaGamma = particle.p() / mass;
    const double rxyz = -betaGamma * ctau * std::log(1 - u);

    if (!charge) {
      decayVertex[0] = particle.vx() + rxyz * (particle.px() / particle.p());
      decayVertex[1] = particle.vy() + rxyz * (particle.py() / particle.p());
      decayVertex[2] = particle.vz() + rxyz * (particle.pz() / particle.p());
    } else {
      o2::math_utils::CircleXYf_t circle;
      o2::track::TrackParCov track = o2::upgrade::convertMCParticleToO2Track(particle, pdgDB);

      float sna{}, csa{};
      track.getCircleParams(mBz, circle, sna, csa);
      const double rxy = rxyz / std::sqrt(1. + track.getTgl() * track.getTgl());
      mTheta = rxy / circle.rC;

      decayVertex[0] = ((particle.vx() - circle.xC) * std::cos(mTheta) - (particle.vy() - circle.yC) * std::sin(mTheta)) + circle.xC;
      decayVertex[1] = ((particle.vy() - circle.yC) * std::cos(mTheta) + (particle.vx() - circle.xC) * std::sin(mTheta)) + circle.yC;
      decayVertex[2] = particle.vz() + rxyz * (particle.pz() / track.getP());
    }
    return decayVertex;
  }

  // Setters
  void setBField(const double b) { mBz = b; }
  void setSeed(const int seed)
  {
    mRand3.SetSeed(seed);   // For decay length sampling
    gRandom->SetSeed(seed); // For TGenPhaseSpace
  }

  // Getters
  [[nodiscard]] float getSecondaryVertexX() const { return static_cast<float>(mVx); }
  [[nodiscard]] float getSecondaryVertexY() const { return static_cast<float>(mVy); }
  [[nodiscard]] float getSecondaryVertexZ() const { return static_cast<float>(mVz); }
  [[nodiscard]] float getDecayRadius() const { return static_cast<float>(std::hypot(mVx, mVy)); }

 private:
  double mBz{20.}; // kG
  double mVx{-1.}, mVy{-1.}, mVz{-1.};
  double mTheta{};
  TRandom3 mRand3;
};

} // namespace o2::upgrade

#endif // ALICE3_CORE_DECAYER_H_
