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
///

#ifndef ALICE3_CORE_DECAYER_H_
#define ALICE3_CORE_DECAYER_H_

#include "ALICE3/Core/TrackUtilities.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TDecayChannel.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TParticlePDG.h>
#include <TRandom3.h>

#include <string>
#include <unordered_map>
#include <vector>
#include <array>
#include <cmath>


namespace o2
{
namespace upgrade
{

class Decayer {
 public:
  // Default constructor
  Decayer() = default;

  template <typename TDatabase>
  std::vector<o2::upgrade::OTFParticle> decayParticle(const TDatabase& pdgDB, const o2::track::TrackParCov& track, const int pdgCode)
  {
    const auto& particleInfo = pdgDB->GetParticle(pdgCode);
    const int charge = particleInfo->Charge() / 3;
    const double mass = particleInfo->Mass();
    std::array<float, 3> mom;
    std::array<float, 3> pos;
    track.getPxPyPzGlo(mom);
    track.getXYZGlo(pos);

    const double u = mRand3.Uniform(0, 1);
    const double ctau = o2::constants::physics::LightSpeedCm2S * particleInfo->Lifetime(); // cm
    const double betaGamma = track.getP() / mass;
    const double rxyz = -betaGamma * ctau * std::log(1 - u);
    double vx, vy, vz;
    double px, py, e;

    if (!charge) {
      vx = pos[0] + rxyz * (mom[0] / track.getP());
      vy = pos[1] + rxyz * (mom[1] / track.getP());
      vz = pos[2] + rxyz * (mom[2] / track.getP());
      px = mom[0];
      py = mom[2];
    } else {
      float sna, csa;
      o2::math_utils::CircleXYf_t circle;
      track.getCircleParams(mBz, circle, sna, csa);
      const double rxy = rxyz / std::sqrt(1. + track.getTgl() * track.getTgl());
      const double theta = rxy / circle.rC;
      
      vx = ((pos[0] - circle.xC) * std::cos(theta) - (pos[1] - circle.yC) * std::sin(theta)) + circle.xC;
      vy = ((pos[1] - circle.yC) * std::cos(theta) + (pos[0] - circle.xC) * std::sin(theta)) + circle.yC;
      vz = mom[2] + rxyz * (mom[2] / track.getP());
      
      px = mom[0] * std::cos(theta) - mom[1] * std::sin(theta);
      py = mom[1] * std::cos(theta) + mom[0] * std::sin(theta);
    }

    double brTotal = 0.;
    e = std::sqrt(mass * mass + px * px + py * py + mom[2] * mom[2]);
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

    TLorentzVector tlv(px, py, mom[2], e);
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
    mRand3.SetSeed(seed); // For decay length sampling
    gRandom->SetSeed(seed); // For TGenPhaseSpace
  }

 void setBField(const double b) { mBz = b; }

 private:
  TRandom3 mRand3;
  double mBz;
};

} // namespace fastsim
} // namespace o2

#endif // ALICE3_CORE_DECAYER_H_
