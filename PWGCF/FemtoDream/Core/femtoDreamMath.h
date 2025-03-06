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

/// \file FemtoDreamMath.h
/// \brief Definition of the FemtoDreamMath Container for math calculations of quantities related to pairs
/// \author Valentina Mantovani Sarti, TU München, valentina.mantovani-sarti@tum.de, Laura Serksnyte, TU München, laura.serksnyte@cern.ch

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMMATH_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMMATH_H_

#include <iostream>

#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TLorentzVector.h"
#include "TMath.h"

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamMath
/// \brief Container for math calculations of quantities related to pairs
class FemtoDreamMath
{
 public:
  /// Compute the k* of a pair of particles
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  template <typename T1, typename T2>
  static float getkstar(const T1& part1, const float mass1, const T2& part2, const float mass2)
  {
    const ROOT::Math::PtEtaPhiMVector vecpart1(part1.pt(), part1.eta(), part1.phi(), mass1);
    const ROOT::Math::PtEtaPhiMVector vecpart2(part2.pt(), part2.eta(), part2.phi(), mass2);
    const ROOT::Math::PtEtaPhiMVector trackSum = vecpart1 + vecpart2;

    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());

    ROOT::Math::PxPyPzMVector PartOneCMS(vecpart1);
    ROOT::Math::PxPyPzMVector PartTwoCMS(vecpart2);

    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);

    const ROOT::Math::PxPyPzMVector trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }
  /// Compute the qij of a pair of particles
  /// \tparam T type of tracks
  /// \param vecparti Particle i PxPyPzMVector
  /// \param vecpartj Particle j PxPyPzMVector
  // The q12 components can be calculated as:
  //             q^mu = (p1-p2)^mu /2 - ((p1-p2)*P/(2P^2))*P^mu
  //             where P = p1+p2
  // Reference: https://www.annualreviews.org/doi/pdf/10.1146/annurev.nucl.55.090704.151533
  // In the following code the above written equation will be expressed as:
  //             q = trackDifference/2 -  scaling * trackSum
  // where scaling is a float number:
  //             scaling = trackDifference*trackSum/(2*trackSum^2) = ((p1-p2)*P/(2P^2))
  // We don't use the reduced vector - no division by 2
  template <typename T>
  static ROOT::Math::PxPyPzEVector getqij(const T& vecparti, const T& vecpartj)
  {
    ROOT::Math::PxPyPzEVector trackSum = vecparti + vecpartj;
    ROOT::Math::PxPyPzEVector trackDifference = vecparti - vecpartj;
    float scaling = trackDifference.Dot(trackSum) / trackSum.Dot(trackSum);
    return trackDifference - scaling * trackSum;
  }

  /// Compute the Q3 of a triplet of particles
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  /// \param part3 Particle 3
  /// \param mass3 Mass of particle 3
  template <typename T>
  static float getQ3(const T& part1, const float mass1, const T& part2, const float mass2, const T& part3, const float mass3)
  {

    const ROOT::Math::PtEtaPhiMVector vecpart01(part1.pt(), part1.eta(), part1.phi(), mass1);
    const ROOT::Math::PtEtaPhiMVector vecpart02(part2.pt(), part2.eta(), part2.phi(), mass2);
    const ROOT::Math::PtEtaPhiMVector vecpart03(part3.pt(), part3.eta(), part3.phi(), mass3);

    const ROOT::Math::PxPyPzEVector vecpart1(vecpart01);
    const ROOT::Math::PxPyPzEVector vecpart2(vecpart02);
    const ROOT::Math::PxPyPzEVector vecpart3(vecpart03);

    ROOT::Math::PxPyPzEVector q12 = getqij(vecpart1, vecpart2);
    ROOT::Math::PxPyPzEVector q23 = getqij(vecpart2, vecpart3);
    ROOT::Math::PxPyPzEVector q31 = getqij(vecpart3, vecpart1);

    float Q32 = q12.M2() + q23.M2() + q31.M2();

    return sqrt(-Q32);
  }

  /// Compute the transverse momentum of a pair of particles
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  template <typename T1, typename T2>
  static float getkT(const T1& part1, const float mass1, const T2& part2, const float mass2)
  {
    const ROOT::Math::PtEtaPhiMVector vecpart1(part1.pt(), part1.eta(), part1.phi(), mass1);
    const ROOT::Math::PtEtaPhiMVector vecpart2(part2.pt(), part2.eta(), part2.phi(), mass2);
    const ROOT::Math::PtEtaPhiMVector trackSum = vecpart1 + vecpart2;
    return 0.5 * trackSum.Pt();
  }

  /// Compute the transverse mass of a pair of particles
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  template <typename T1, typename T2>
  static float getmT(const T1& part1, const float mass1, const T2& part2, const float mass2)
  {
    return std::sqrt(std::pow(getkT(part1, mass1, part2, mass2), 2.) + std::pow(0.5 * (mass1 + mass2), 2.));
  }

  template <typename T1>
  static float getInvMassCascade(const T1& trackpos, const float masspos, const T1& trackneg, const float massneg, const T1& trackbach, const float massbach, const float massv0)
  {
    // calculate the invariant mass
    const ROOT::Math::PtEtaPhiMVector posDaug(trackpos.pt(), trackpos.eta(), trackpos.phi(), masspos);
    const ROOT::Math::PtEtaPhiMVector negDaug(trackneg.pt(), trackneg.eta(), trackneg.phi(), massneg);
    const ROOT::Math::PtEtaPhiMVector bachDaug(trackbach.pt(), trackbach.eta(), trackbach.phi(), massbach);
    const ROOT::Math::PxPyPzMVector v0(posDaug.Px() + negDaug.Px(), posDaug.Py() + negDaug.Py(), posDaug.Pz() + negDaug.Pz(), massv0);
    const ROOT::Math::PxPyPzMVector casc = v0 + bachDaug;

    return casc.M();
  }
};

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMMATH_H_
