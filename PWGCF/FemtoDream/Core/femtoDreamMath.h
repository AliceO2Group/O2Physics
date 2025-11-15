// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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

#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include <Math/VectorUtil.h>

#include <iostream>
#include <vector>

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

  /// Compute the 3d components of the pair momentum in LCMS and PRF
  /// Copy from femto universe
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  /// \param isiden Identical or non-identical particle pair
  template <typename T>
  static std::vector<double> newpairfunc(const T& part1, const float mass1, const T& part2, const float mass2, bool isiden)
  {
    const double e1 = std::sqrt(std::pow(part1.px(), 2) + std::pow(part1.py(), 2) + std::pow(part1.pz(), 2) + std::pow(mass1, 2));
    const double e2 = std::sqrt(std::pow(part2.px(), 2) + std::pow(part2.py(), 2) + std::pow(part2.pz(), 2) + std::pow(mass2, 2));

    const ROOT::Math::PxPyPzEVector vecpart1(part1.px(), part1.py(), part1.pz(), e1);
    const ROOT::Math::PxPyPzEVector vecpart2(part2.px(), part2.py(), part2.pz(), e2);
    const ROOT::Math::PxPyPzEVector trackSum = vecpart1 + vecpart2;

    std::vector<double> vect;

    const double tPx = trackSum.px();
    const double tPy = trackSum.py();
    const double tPz = trackSum.pz();
    const double tE = trackSum.E();

    const double tPtSq = (tPx * tPx + tPy * tPy);
    const double tMtSq = (tE * tE - tPz * tPz);
    const double tM = std::sqrt(tMtSq - tPtSq);
    const double tMt = std::sqrt(tMtSq);
    const double tPt = std::sqrt(tPtSq);

    // Boost to LCMS

    const double beta_LCMS = tPz / tE;
    const double gamma_LCMS = tE / tMt;

    const double fDKOut = (part1.px() * tPx + part1.py() * tPy) / tPt;
    const double fDKSide = (-part1.px() * tPy + part1.py() * tPx) / tPt;
    const double fDKLong = gamma_LCMS * (part1.pz() - beta_LCMS * e1);
    const double fDE = gamma_LCMS * (e1 - beta_LCMS * part1.pz());

    const double px1LCMS = fDKOut;
    const double py1LCMS = fDKSide;
    const double pz1LCMS = fDKLong;
    const double pE1LCMS = fDE;

    const double px2LCMS = (part2.px() * tPx + part2.py() * tPy) / tPt;
    const double py2LCMS = (part2.py() * tPx - part2.px() * tPy) / tPt;
    const double pz2LCMS = gamma_LCMS * (part2.pz() - beta_LCMS * e2);
    const double pE2LCMS = gamma_LCMS * (e2 - beta_LCMS * part2.pz());

    const double fDKOutLCMS = px1LCMS - px2LCMS;
    const double fDKSideLCMS = py1LCMS - py2LCMS;
    const double fDKLongLCMS = pz1LCMS - pz2LCMS;

    // Boost to PRF

    const double betaOut = tPt / tMt;
    const double gammaOut = tMt / tM;

    const double fDKOutPRF = gammaOut * (fDKOutLCMS - betaOut * (pE1LCMS - pE2LCMS));
    const double fDKSidePRF = fDKSideLCMS;
    const double fDKLongPRF = fDKLongLCMS;
    const double fKOut = gammaOut * (fDKOut - betaOut * fDE);

    const double qlcms = std::sqrt(fDKOutLCMS * fDKOutLCMS + fDKSideLCMS * fDKSideLCMS + fDKLongLCMS * fDKLongLCMS);
    const double qinv = std::sqrt(fDKOutPRF * fDKOutPRF + fDKSidePRF * fDKSidePRF + fDKLongPRF * fDKLongPRF);
    const double kstar = std::sqrt(fKOut * fKOut + fDKSide * fDKSide + fDKLong * fDKLong);

    if (isiden) {
      vect.push_back(qinv);
      vect.push_back(fDKOutLCMS);
      vect.push_back(fDKSideLCMS);
      vect.push_back(fDKLongLCMS);
      vect.push_back(qlcms);
    } else {
      vect.push_back(kstar);
      vect.push_back(fDKOut);
      vect.push_back(fDKSide);
      vect.push_back(fDKLong);
    }
    return vect;
  }

  /// Compute the phi angular of a pair with respect to the event plane
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  template <typename T1, typename T2>
  static float getPairPhiEP(const T1& part1, const float mass1, const T2& part2, const float mass2, const float Psi_ep)
  {
    const ROOT::Math::PtEtaPhiMVector vecpart1(part1.pt(), part1.eta(), part1.phi(), mass1);
    const ROOT::Math::PtEtaPhiMVector vecpart2(part2.pt(), part2.eta(), part2.phi(), mass2);
    const ROOT::Math::PtEtaPhiMVector trackSum = vecpart1 + vecpart2;
    float phi_pair_onPsi = TVector2::Phi_mpi_pi(trackSum.Phi() - Psi_ep);
    phi_pair_onPsi = TMath::Abs(phi_pair_onPsi);
    return phi_pair_onPsi;
  }

  /// Compute the phi angular of a pair according to plane-calibarated second particle
  template <typename T1, typename T2>
  static float getPairPhiEP(const T1& part1, const float mass1, const T2& part2, const float mass2, const float Psi_ep1, const float Psi_ep2)
  {
    const ROOT::Math::PtEtaPhiMVector vecpart1(part1.pt(), part1.eta(), part1.phi(), mass1);
    const ROOT::Math::PtEtaPhiMVector vecpart2(part2.pt(), part2.eta(), part2.phi(), mass2);
    const float psidiff = Psi_ep2 - Psi_ep1;
    // rotate phi of part2
    const float newPhi2 = TVector2::Phi_mpi_pi(vecpart2.Phi() - psidiff);
    const ROOT::Math::PtEtaPhiMVector vecpart2_calibd(vecpart2.Pt(), vecpart2.Eta(), newPhi2, vecpart2.M());
    const ROOT::Math::PtEtaPhiMVector trackSum = vecpart1 + vecpart2_calibd;
    // reCalibrate part2 phi with respect to part1 event plane
    float phi_pair_onPsi = TVector2::Phi_mpi_pi(trackSum.Phi() - Psi_ep1);
    phi_pair_onPsi = TMath::Abs(phi_pair_onPsi);
    return phi_pair_onPsi;
  }
};

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMMATH_H_
