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

/// \brief Functions which cut on particle pairs (decays, conversions, two-track cuts) adapted for data from UD tables
/// Based on the code "PWGCF/Core/PairCuts.h" made by Jan Fiete Grosse-Oetringhaus
/// Author:

#ifndef PWGUD_CORE_UPCPAIRCUTS_H_
#define PWGUD_CORE_UPCPAIRCUTS_H_

#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include <cmath>

class UPCPairCuts
{
 public:
  enum Particle { Photon = 0,
                  K0,
                  Lambda,
                  Phi,
                  Rho,
                  ParticlesLastEntry };

  void setHistogramRegistry(o2::framework::HistogramRegistry* registry) { histogramRegistry = registry; }

  void setPairCut(Particle particle, float cut)
  {
    LOGF(info, "Enabled pair cut for %d with value %f", static_cast<int>(particle), cut);
    mCuts[particle] = cut;
    if (histogramRegistry != nullptr && histogramRegistry->contains(HIST("ControlConvResonances")) == false) {
      histogramRegistry->add("ControlConvResonances", "", {o2::framework::HistType::kTH2F, {{6, -0.5, 5.5, "id"}, {500, -0.5, 0.5, "delta mass"}}});
    }
  }

  void setTwoTrackCuts(float distance = 0.02f, float radius = 0.8f)
  {
    LOGF(info, "Enabled two-track cut with distance %f and radius %f", distance, radius);
    mTwoTrackDistance = distance;
    mTwoTrackRadius = radius;

    if (histogramRegistry != nullptr && histogramRegistry->contains(HIST("TwoTrackDistancePt_0")) == false) {
      histogramRegistry->add("TwoTrackDistancePt_0", "", {o2::framework::HistType::kTH3F, {{100, -0.15, 0.15, "#Delta#eta"}, {100, -0.05, 0.05, "#Delta#varphi^{*}_{min}"}, {20, 0, 10, "#Delta p_{T}"}}});
      histogramRegistry->addClone("TwoTrackDistancePt_0", "TwoTrackDistancePt_1");
    }
  }

  template <typename T>
  bool conversionCuts(T const& track1, T const& track2);

  template <typename T>
  bool twoTrackCut(T const& track1, T const& track2);

 protected:
  float mCuts[ParticlesLastEntry] = {-1};
  float mTwoTrackDistance = -1; // distance below which the pair is flagged as to be removed
  float mTwoTrackRadius = 0.8f; // radius at which the two track cuts are applied
  int magField = 5;             // magField: B field in kG

  o2::framework::HistogramRegistry* histogramRegistry = nullptr; // if set, control histograms are stored here

  template <typename T>
  bool conversionCut(T const& track1, T const& track2, Particle conv, double cut);

  template <typename T>
  double getInvMassSquared(T const& track1, double m0_1, T const& track2, double m0_2);

  template <typename T>
  double getInvMassSquaredFast(T const& track1, double m0_1, T const& track2, double m0_2);

  template <typename T>
  float getDPhiStar(T const& track1, T const& track2, float radius, int magField);
};

template <typename T>
bool UPCPairCuts::conversionCuts(T const& track1, T const& track2)
{
  // skip if like sign
  if (track1.sign() * track2.sign() > 0) {
    return false;
  }

  for (int i = 0; i < static_cast<int>(ParticlesLastEntry); i++) {
    Particle particle = static_cast<Particle>(i);
    if (mCuts[i] > 0) {
      if (conversionCut(track1, track2, particle, mCuts[i])) {
        return true;
      }
      if (particle == Lambda) {
        if (conversionCut(track2, track1, particle, mCuts[i])) {
          return true;
        }
      }
    }
  }

  return false;
}

template <typename T>
bool UPCPairCuts::twoTrackCut(T const& track1, T const& track2)
{
  // the variables & cut have been developed in Run 1 by the CF - HBT group
  //
  // Parameters:
  //   magField: B field in kG

  auto deta = eta(track1.px(), track1.py(), track1.pz()) - eta(track2.px(), track2.py(), track2.pz());

  // optimization
  if (std::fabs(deta) < mTwoTrackDistance * 2.5 * 3) {
    // check first boundaries to see if is worth to loop and find the minimum
    float dphistar1 = getDPhiStar(track1, track2, mTwoTrackRadius, magField);
    float dphistar2 = getDPhiStar(track1, track2, 2.5, magField);

    const float kLimit = mTwoTrackDistance * 3;

    if (std::fabs(dphistar1) < kLimit || std::fabs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0) {
      float dphistarminabs = 1e5;
      float dphistarmin = 1e5;
      for (Double_t rad = mTwoTrackRadius; rad < 2.51; rad += 0.01) {
        float dphistar = getDPhiStar(track1, track2, rad, magField);

        float dphistarabs = std::fabs(dphistar);

        if (dphistarabs < dphistarminabs) {
          dphistarmin = dphistar;
          dphistarminabs = dphistarabs;
        }
      }

      if (histogramRegistry != nullptr) {
        histogramRegistry->fill(HIST("TwoTrackDistancePt_0"), deta, dphistarmin, std::fabs(track1.pt() - track2.pt()));
      }

      if (dphistarminabs < mTwoTrackDistance && std::fabs(deta) < mTwoTrackDistance) {
        // LOGF(debug, "Removed track pair %ld %ld with %f %f %f %f %d %f %f %d %d", track1.index(), track2.index(), deta, dphistarminabs, track1.phi2(), track1.pt(), track1.sign(), track2.phi2(), track2.pt(), track2.sign(), magField);
        return true;
      }

      if (histogramRegistry != nullptr) {
        histogramRegistry->fill(HIST("TwoTrackDistancePt_1"), deta, dphistarmin, std::fabs(track1.pt() - track2.pt()));
      }
    }
  }

  return false;
}

template <typename T>
bool UPCPairCuts::conversionCut(T const& track1, T const& track2, Particle conv, double cut)
{
  // LOGF(info, "pt is %f %f", track1.pt(), track2.pt());

  if (cut < 0) {
    return false;
  }

  double massD1, massD2, massM;

  switch (conv) {
    case Photon:
      massD1 = o2::constants::physics::MassElectron;
      massD2 = o2::constants::physics::MassElectron;
      massM = 0;
      break;
    case K0:
      massD1 = o2::constants::physics::MassPiPlus;
      massD2 = o2::constants::physics::MassPiPlus;
      massM = o2::constants::physics::MassK0;
      break;
    case Lambda:
      massD1 = o2::constants::physics::MassProton;
      massD2 = o2::constants::physics::MassPiPlus;
      massM = o2::constants::physics::MassLambda0;
      break;
    case Phi:
      massD1 = o2::constants::physics::MassKPlus;
      massD2 = o2::constants::physics::MassKPlus;
      massM = o2::constants::physics::MassPhi;
      break;
    case Rho:
      massD1 = o2::constants::physics::MassPiPlus;
      massD2 = o2::constants::physics::MassPiPlus;
      massM = 0.770;
      break;
    default:
      LOGF(fatal, "Particle now known");
      return false;
      break;
  }

  auto massC = getInvMassSquaredFast(track1, massD1, track2, massD2);

  if (std::fabs(massC - massM * massM) > cut * 5) {
    return false;
  }

  massC = getInvMassSquared(track1, massD1, track2, massD2);

  if (histogramRegistry != nullptr) {
    histogramRegistry->fill(HIST("ControlConvResonances"), static_cast<int>(conv), massC - massM * massM);
  }

  if (massC > (massM - cut) * (massM - cut) && massC < (massM + cut) * (massM + cut)) {
    return true;
  }

  return false;
}

template <typename T>
double UPCPairCuts::getInvMassSquared(T const& track1, double m0_1, T const& track2, double m0_2)
{
  // calculate inv mass squared
  // same can be achieved, but with more computing time with
  /*TLorentzVector photon, p1, p2;
  p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
  p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
  photon = p1+p2;
  photon.M()*/

  float tantheta1 = 1e10;

  if (eta(track1.px(), track1.py(), track1.pz()) < -1e-10 || eta(track1.px(), track1.py(), track1.pz()) > 1e-10) {
    float expTmp = std::exp(-eta(track1.px(), track1.py(), track1.pz()));
    tantheta1 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
  }

  float tantheta2 = 1e10;
  if (eta(track2.px(), track2.py(), track2.pz()) < -1e-10 || eta(track2.px(), track2.py(), track2.pz()) > 1e-10) {
    float expTmp = std::exp(-eta(track2.px(), track2.py(), track2.pz()));
    tantheta2 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
  }

  float e1squ = m0_1 * m0_1 + track1.pt() * track1.pt() * (1.0 + 1.0 / tantheta1 / tantheta1);
  float e2squ = m0_2 * m0_2 + track2.pt() * track2.pt() * (1.0 + 1.0 / tantheta2 / tantheta2);

  float mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * (std::sqrt(e1squ * e2squ) - (track1.pt() * track2.pt() * (std::cos(phi(track1.px(), track1.py()) - phi(track2.px(), track2.py())) + 1.0 / tantheta1 / tantheta2)));

  // LOGF(debug, "%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2);

  return mass2;
}

template <typename T>
double UPCPairCuts::getInvMassSquaredFast(T const& track1, double m0_1, T const& track2, double m0_2)
{
  // calculate inv mass squared approximately

  const float eta1 = eta(track1.px(), track1.py(), track1.pz());
  const float eta2 = eta(track2.px(), track2.py(), track2.pz());
  const float phi1 = phi(track1.px(), track1.py());
  const float phi2 = phi(track2.px(), track2.py());
  const float pt1 = track1.pt();
  const float pt2 = track2.pt();

  float tantheta1 = 1e10f;

  if (eta1 < -1e-10f || eta1 > 1e-10f) {
    float expTmp = 1.0f - eta1 + eta1 * eta1 / 2.0f - eta1 * eta1 * eta1 / 6.0f + eta1 * eta1 * eta1 * eta1 / 24.0f;
    tantheta1 = 2.0f * expTmp / (1.0f - expTmp * expTmp);
  }

  float tantheta2 = 1e10f;
  if (eta2 < -1e-10f || eta2 > 1e-10f) {
    float expTmp = 1.0f - eta2 + eta2 * eta2 / 2.0f - eta2 * eta2 * eta2 / 6.0f + eta2 * eta2 * eta2 * eta2 / 24.0f;
    tantheta2 = 2.0f * expTmp / (1.0f - expTmp * expTmp);
  }

  float e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0f + 1.0f / tantheta1 / tantheta1);
  float e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0f + 1.0f / tantheta2 / tantheta2);

  // fold onto 0...pi
  float deltaPhi = std::fabs(phi1 - phi2);
  while (deltaPhi > o2::constants::math::TwoPI) {
    deltaPhi -= o2::constants::math::TwoPI;
  }
  if (deltaPhi > o2::constants::math::PI) {
    deltaPhi = o2::constants::math::TwoPI - deltaPhi;
  }

  float cosDeltaPhi = 0;
  if (deltaPhi < o2::constants::math::PI / 3.0f) {
    cosDeltaPhi = 1.0 - deltaPhi * deltaPhi / 2 + deltaPhi * deltaPhi * deltaPhi * deltaPhi / 24;
  } else if (deltaPhi < 2.0f * o2::constants::math::PI / 3.0f) {
    cosDeltaPhi = -(deltaPhi - o2::constants::math::PI / 2) + 1.0 / 6 * std::pow((deltaPhi - o2::constants::math::PI / 2), 3);
  } else {
    cosDeltaPhi = -1.0f + 1.0f / 2.0f * (deltaPhi - o2::constants::math::PI) * (deltaPhi - o2::constants::math::PI) - 1.0f / 24.0f * std::pow(deltaPhi - o2::constants::math::PI, 4.0f);
  }

  double mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2.0f * (std::sqrt(e1squ * e2squ) - (pt1 * pt2 * (cosDeltaPhi + 1.0f / tantheta1 / tantheta2)));

  // LOGF(debug, "%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2);

  return mass2;
}

template <typename T>
float UPCPairCuts::getDPhiStar(T const& track1, T const& track2, float radius, int magField)
{
  //
  // calculates dphistar
  //

  auto phi1 = phi(track1.px(), track1.py());
  auto pt1 = track1.pt();
  auto charge1 = track1.sign();

  auto phi2 = phi(track2.px(), track2.py());
  auto pt2 = track2.pt();
  auto charge2 = track2.sign();

  float dphistar = phi1 - phi2 - charge1 * std::asin(0.015 * magField * radius / pt1) + charge2 * std::asin(0.015 * magField * radius / pt2);

  if (dphistar > o2::constants::math::PI) {
    dphistar = o2::constants::math::TwoPI - dphistar;
  }
  if (dphistar < -o2::constants::math::PI) {
    dphistar = -o2::constants::math::TwoPI - dphistar;
  }
  if (dphistar > o2::constants::math::PI) { // might look funny but is needed
    dphistar = o2::constants::math::TwoPI - dphistar;
  }

  return dphistar;
}

#endif // PWGUD_CORE_UPCPAIRCUTS_H_
