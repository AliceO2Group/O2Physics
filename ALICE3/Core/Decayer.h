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

#include "ALICE3/Core/TrackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/Track.h>

#include <TDecayChannel.h> // IWYU pragma: keep
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace o2
{
namespace upgrade
{

class Decayer
{
 public:
  // Default constructor
  Decayer() = default;

  template <typename TDatabase>
  std::vector<o2::upgrade::OTFParticle> decayParticle(const TDatabase& pdgDB, const OTFParticle& particle)
  {
    const auto& particleInfo = pdgDB->GetParticle(particle.pdgCode());
    if (!particleInfo) {
      return {};
    }

    const int charge = particleInfo->Charge() / 3;
    const double mass = particleInfo->Mass();

    const double u = mRand3.Uniform(0.001, 0.999);
    const double ctau = o2::constants::physics::LightSpeedCm2S * particleInfo->Lifetime(); // cm
    const double betaGamma = particle.p() / mass;
    const double rxyz = -betaGamma * ctau * std::log(1 - u);
    double vx, vy, vz;
    double px, py, e;

    if (!charge) {
      vx = particle.vx() + rxyz * (particle.px() / particle.p());
      vy = particle.vy() + rxyz * (particle.py() / particle.p());
      vz = particle.vz() + rxyz * (particle.pz() / particle.p());
      px = particle.px();
      py = particle.py();
    } else {
      o2::track::TrackParCov track;
      o2::math_utils::CircleXYf_t circle;
      o2::upgrade::convertOTFParticleToO2Track(particle, track, pdgDB);

      float sna, csa;
      track.getCircleParams(mBz, circle, sna, csa);
      const double rxy = rxyz / std::sqrt(1. + track.getTgl() * track.getTgl());
      const double theta = rxy / circle.rC;

      vx = ((particle.vx() - circle.xC) * std::cos(theta) - (particle.vy() - circle.yC) * std::sin(theta)) + circle.xC;
      vy = ((particle.vy() - circle.yC) * std::cos(theta) + (particle.vx() - circle.xC) * std::sin(theta)) + circle.yC;
      vz = particle.vz() + rxyz * (particle.pz() / track.getP());

      px = particle.px() * std::cos(theta) - particle.py() * std::sin(theta);
      py = particle.py() * std::cos(theta) + particle.px() * std::sin(theta);
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
      particle.setVxVyVz(vx, vy, vz);
      particle.setPxPyPzE(dau.Px(), dau.Py(), dau.Pz(), dau.E());
      decayProducts.push_back(particle);
    }

    return decayProducts;
  }

  void setSeed(const int seed)
  {
    mRand3.SetSeed(seed);   // For decay length sampling
    gRandom->SetSeed(seed); // For TGenPhaseSpace
  }

  void setBField(const double b) { mBz = b; }

 private:
  TRandom3 mRand3;
  double mBz;
};

} // namespace upgrade
} // namespace o2

#endif // ALICE3_CORE_DECAYER_H_
