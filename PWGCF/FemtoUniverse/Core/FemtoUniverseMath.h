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

/// \file FemtoUniverseMath.h
/// \brief Definition of the FemtoUniverseMath Container for math calculations of quantities related to pairs
/// \author Valentina Mantovani Sarti, TU München, valentina.mantovani-sarti@tum.de
/// \author Laura Serksnyte, TU München, laura.serksnyte@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEMATH_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEMATH_H_

#include <iostream>
#include <vector>
#include <algorithm>

#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TLorentzVector.h"
#include "TMath.h"

namespace o2::analysis::femtoUniverse
{

/// \class FemtoUniverseMath
/// \brief Container for math calculations of quantities related to pairs
class FemtoUniverseMath
{
 public:
  /// Compute the k* of a pair of particles
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  template <typename T>
  static float getkstar(const T& part1, const float mass1, const T& part2, const float mass2)
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
    float E1 = sqrt(pow(part1.px(), 2) + pow(part1.py(), 2) + pow(part1.pz(), 2) + pow(mass1, 2));
    float E2 = sqrt(pow(part2.px(), 2) + pow(part2.py(), 2) + pow(part2.pz(), 2) + pow(mass2, 2));
    float E3 = sqrt(pow(part3.px(), 2) + pow(part3.py(), 2) + pow(part3.pz(), 2) + pow(mass3, 2));

    const ROOT::Math::PxPyPzEVector vecpart1(part1.px(), part1.py(), part1.pz(), E1);
    const ROOT::Math::PxPyPzEVector vecpart2(part2.px(), part2.py(), part2.pz(), E2);
    const ROOT::Math::PxPyPzEVector vecpart3(part3.px(), part3.py(), part3.pz(), E3);

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
  template <typename T>
  static float getkT(const T& part1, const float mass1, const T& part2, const float mass2)
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
  template <typename T>
  static float getmT(const T& part1, const float mass1, const T& part2, const float mass2)
  {
    return std::sqrt(std::pow(getkT(part1, mass1, part2, mass2), 2.) + std::pow(0.5 * (mass1 + mass2), 2.));
  }

  /// Compute the 3d components of the pair momentum in LCMS and PRF
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  /// \param isiden Identical or non-identical particle pair
  /// \param islcms LCMS or PRF
  template <typename T>
  static std::vector<double> getpairmom3d(const T& part1, const float mass1, const T& part2, const float mass2, bool isiden, bool islcms)
  {
    const double E1 = sqrt(pow(part1.px(), 2) + pow(part1.py(), 2) + pow(part1.pz(), 2) + pow(mass1, 2));
    const double E2 = sqrt(pow(part2.px(), 2) + pow(part2.py(), 2) + pow(part2.pz(), 2) + pow(mass2, 2));

    const ROOT::Math::PxPyPzEVector vecpart1(part1.px(), part1.py(), part1.pz(), E1);
    const ROOT::Math::PxPyPzEVector vecpart2(part2.px(), part2.py(), part2.pz(), E2);
    const ROOT::Math::PxPyPzEVector trackSum = vecpart1 + vecpart2;

    std::vector<double> vect;

    const double tPx = trackSum.px();
    const double tPy = trackSum.py();
    const double tPz = trackSum.pz();
    const double tPE = trackSum.E();

    const double tPt = trackSum.pt();
    const double tMt = trackSum.mt();
    const double tPinv = std::sqrt((tMt * tMt) - (tPt * tPt));

    float nullmass = 0.0;
    const double m1 = std::max(nullmass, mass1);
    const double m2 = std::max(nullmass, mass2);

    const double tQinvL = std::pow((E1 - E2), 2) - std::pow((part1.px() - part2.px()), 2) -
                          std::pow((part1.py() - part2.py()), 2) - std::pow((part1.pz() - part2.pz()), 2);

    double tQ = (m1 - m2) / tPinv;
    tQ = ::sqrt(tQ * tQ - tQinvL);

    const double fKStarCalc = tQ / 2.0;
    vect.push_back(fKStarCalc);

    // Boost to LCMS

    const double beta = tPz / tPE;
    const double gamma = tPE / tMt;

    const double px1L = (part1.px() * tPx + part1.py() * tPy) / tPt;
    const double py1L = (-part1.px() * tPy + part1.py() * tPx) / tPt;
    const double pz1L = gamma * (part1.pz() - beta * E1);
    const double pE1L = gamma * (E1 - beta * part1.pz());

    const double px2L = (part2.px() * tPx + part2.py() * tPy) / tPt;
    const double py2L = (-part2.px() * tPy + part2.py() * tPx) / tPt;
    const double pz2L = gamma * (part2.pz() - beta * E2);
    const double pE2L = gamma * (E2 - beta * part2.pz());

    double fDKOutLCMS;
    double fDKSideLCMS;
    double fDKLongLCMS;

    double fDKOutPRF;
    double fDKSidePRF;
    double fDKLongPRF;

    if (!isiden) {
      fDKOutLCMS = px1L;
      fDKSideLCMS = py1L;
      fDKLongLCMS = pz1L;
    } else {
      fDKOutLCMS = px1L - px2L;
      fDKSideLCMS = py1L - py2L;
      fDKLongLCMS = pz1L - pz2L;
    }

    // Boost to PRF
    const double betaOut = tPt / tMt;
    const double gammaOut = tMt / tPinv;

    if (!isiden) {
      fDKOutPRF = gammaOut * (fDKOutLCMS - betaOut * pE1L);
      fDKSidePRF = fDKSideLCMS;
      fDKLongPRF = fDKLongLCMS;
    } else {
      fDKOutPRF = gammaOut * (fDKOutLCMS - betaOut * (pE1L - pE2L));
      fDKSidePRF = fDKSideLCMS;
      fDKLongPRF = fDKLongLCMS;
    }

    if (islcms) {
      vect.push_back(fDKOutLCMS);
      vect.push_back(fDKSideLCMS);
      vect.push_back(fDKLongLCMS);
    } else {
      vect.push_back(fDKOutPRF);
      vect.push_back(fDKSidePRF);
      vect.push_back(fDKLongPRF);
    }
    return vect;
  }
};

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEMATH_H_
