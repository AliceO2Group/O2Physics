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

/// \file RecoDecay.h
/// \brief Implementation of the RecoDecay class
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#ifndef COMMON_CORE_RECODECAY_H_
#define COMMON_CORE_RECODECAY_H_

#include <algorithm> // std::find
#include <array>     // std::array
#include <cmath>     // std::abs, std::sqrt
#include <utility>   // std::move
#include <vector>    // std::vector

#include "CommonConstants/MathConstants.h"

/// Base class for calculating properties of reconstructed decays
///
/// Provides static helper functions for:
/// - useful arithmetic operations and basic vector algebra
/// - calculation of kinematic quantities
/// - calculation of topological properties of secondary vertices
/// - Monte Carlo matching of decays at track and particle level

class RecoDecay
{
 public:
  /// Default constructor
  RecoDecay() = default;

  /// Default destructor
  ~RecoDecay() = default;

  // mapping of charm-hadron origin type
  enum OriginType { None = 0,
                    Prompt,
                    NonPrompt };

  static constexpr int8_t PdgStatusCodeAfterFlavourOscillation = 92; // decay products after B0(s) flavour oscillation

  // Auxiliary functions

  /// Sums numbers.
  /// \param args  arbitrary number of numbers of arbitrary types
  /// \return sum of numbers
  template <typename... T>
  static auto sum(const T&... args)
  {
    return (args + ...);
  }

  /// Squares a number.
  /// \note Promotes number to double before squaring to avoid precision loss in float multiplication.
  /// \param num  a number of arbitrary type
  /// \return number squared
  template <typename T>
  static double sq(T num)
  {
    return static_cast<double>(num) * static_cast<double>(num);
  }

  /// Sums squares of numbers.
  /// \note Promotes numbers to double before squaring to avoid precision loss in float multiplication.
  /// \param args  arbitrary number of numbers of arbitrary types
  /// \return sum of squares of numbers
  template <typename... T>
  static double sumOfSquares(const T&... args)
  {
    return ((static_cast<double>(args) * static_cast<double>(args)) + ...);
  }

  /// Calculates square root of a sum of squares of numbers.
  /// \param args  arbitrary number of numbers of arbitrary types
  /// \return square root of sum of squares of numbers
  template <typename... T>
  static double sqrtSumOfSquares(const T&... args)
  {
    return std::sqrt(sumOfSquares(args...));
  }

  /// Sums i-th elements of containers.
  /// \param index  element index
  /// \param args  pack of containers of elements accessible by index
  /// \return sum of i-th elements
  template <typename... T>
  static auto getElement(int index, const T&... args)
  {
    return (args[index] + ...);
  }

  /// Sums 3-vectors.
  /// \param args  pack of 3-vector arrays
  /// \return sum of vectors
  template <typename... T>
  static auto sumOfVec(const std::array<T, 3>&... args)
  {
    return std::array{getElement(0, args...), getElement(1, args...), getElement(2, args...)};
  }

  /// Calculates scalar product of vectors.
  /// \note Promotes numbers to double to avoid precision loss in float multiplication.
  /// \param N  dimension
  /// \param vec1,vec2  vectors
  /// \return scalar product
  template <std::size_t N, typename T, typename U>
  static double dotProd(const std::array<T, N>& vec1, const std::array<U, N>& vec2)
  {
    double res{0};
    for (std::size_t iDim = 0; iDim < N; ++iDim) {
      res += static_cast<double>(vec1[iDim]) * static_cast<double>(vec2[iDim]);
    }
    return res;
  }

  /// FIXME: probably cross and dot products should be in some utility class
  /// Calculates cross product of vectors in three dimensions.
  /// \note Promotes numbers to double to avoid precision loss in float multiplication.
  /// \param vec1,vec2  vectors
  /// \return cross-product vector
  template <typename T, typename U>
  static std::array<double, 3> crossProd(const std::array<T, 3>& vec1, const std::array<U, 3>& vec2)
  {
    return std::array<double, 3>{(static_cast<double>(vec1[1]) * static_cast<double>(vec2[2])) - (static_cast<double>(vec1[2]) * static_cast<double>(vec2[1])),
                                 (static_cast<double>(vec1[2]) * static_cast<double>(vec2[0])) - (static_cast<double>(vec1[0]) * static_cast<double>(vec2[2])),
                                 (static_cast<double>(vec1[0]) * static_cast<double>(vec2[1])) - (static_cast<double>(vec1[1]) * static_cast<double>(vec2[0]))};
  }

  /// Calculates magnitude squared of a vector.
  /// \param N  dimension
  /// \param vec  vector
  /// \return magnitude squared
  template <std::size_t N, typename T>
  static double mag2(const std::array<T, N>& vec)
  {
    return dotProd(vec, vec);
  }

  /// Calculates 3D distance between two points.
  /// \param point1,point2  {x, y, z} coordinates of points
  /// \return 3D distance between two points
  template <typename T, typename U>
  static double distance(const T& point1, const U& point2)
  {
    return sqrtSumOfSquares(point1[0] - point2[0], point1[1] - point2[1], point1[2] - point2[2]);
  }

  /// Calculates 2D {x, y} distance between two points.
  /// \param point1,point2  {x, y, z} or {x, y} coordinates of points
  /// \return 2D {x, y} distance between two points
  template <typename T, typename U>
  static double distanceXY(const T& point1, const U& point2)
  {
    return sqrtSumOfSquares(point1[0] - point2[0], point1[1] - point2[1]);
  }

  // Calculation of kinematic quantities

  /// Calculates pseudorapidity.
  /// \param mom  3-momentum array
  /// \return pseudorapidity
  template <typename T>
  static double eta(const std::array<T, 3>& mom)
  {
    // eta = arctanh(pz/p)
    if (std::abs(mom[0]) < o2::constants::math::Almost0 && std::abs(mom[1]) < o2::constants::math::Almost0) { // very small px and py
      return static_cast<double>(mom[2] > 0 ? o2::constants::math::VeryBig : -o2::constants::math::VeryBig);
    }
    return static_cast<double>(std::atanh(mom[2] / p(mom)));
  }

  /// Calculates rapidity.
  /// \param mom  3-momentum array
  /// \param mass  mass
  /// \return rapidity
  template <typename T, typename U>
  static double y(const std::array<T, 3>& mom, U mass)
  {
    // y = arctanh(pz/E)
    return std::atanh(mom[2] / e(mom, mass));
  }

  /// Calculates azimuth from x and y components.
  /// \param x,y  {x, y} components
  /// \return azimuth within [0, 2π]
  template <typename T, typename U>
  static double phi(T x, U y)
  {
    // conversion from [-π, +π] returned by atan2 to [0, 2π]
    return std::atan2(static_cast<double>(-y), static_cast<double>(-x)) + o2::constants::math::PI;
  }

  /// Calculates azimuth of a vector.
  /// \note Elements 0 and 1 are expected to represent the x and y vector components, respectively.
  /// \param vec  vector (container of elements accessible by index)
  /// \return azimuth within [0, 2π]
  template <typename T>
  static double phi(const T& vec)
  {
    return phi(vec[0], vec[1]);
  }

  /// Constrains angle to be within a range.
  /// \note Inspired by TVector2::Phi_0_2pi in ROOT.
  /// \param angle  angle
  /// \param min  minimum of the range
  /// \return value within [min, min + 2π).
  template <typename T, typename U = float>
  static T constrainAngle(T angle, U min = 0.)
  {
    while (angle < min) {
      angle += o2::constants::math::TwoPI;
    }
    while (angle >= min + o2::constants::math::TwoPI) {
      angle -= o2::constants::math::TwoPI;
    }
    return (T)angle;
  }

  /// Calculates cosine of pointing angle.
  /// \param posPV  {x, y, z} position of the primary vertex
  /// \param posSV  {x, y, z} position of the secondary vertex
  /// \param mom  3-momentum array
  /// \return cosine of pointing angle
  template <typename T, typename U, typename V>
  static double cpa(const T& posPV, const U& posSV, const std::array<V, 3>& mom)
  {
    // CPA = (l . p)/(|l| |p|)
    auto lineDecay = std::array{posSV[0] - posPV[0], posSV[1] - posPV[1], posSV[2] - posPV[2]};
    auto cos = dotProd(lineDecay, mom) / std::sqrt(mag2(lineDecay) * mag2(mom));
    if (cos < -1.) {
      return -1.;
    }
    if (cos > 1.) {
      return 1.;
    }
    return cos;
  }

  /// Calculates cosine of pointing angle in the {x, y} plane.
  /// \param posPV  {x, y, z} or {x, y} position of the primary vertex
  /// \param posSV  {x, y, z} or {x, y} position of the secondary vertex
  /// \param mom  {x, y, z} or {x, y} momentum array
  /// \return cosine of pointing angle in {x, y}
  template <std::size_t N, typename T, typename U, typename V>
  static double cpaXY(const T& posPV, const U& posSV, const std::array<V, N>& mom)
  {
    // CPAXY = (r . pT)/(|r| |pT|)
    auto lineDecay = std::array{posSV[0] - posPV[0], posSV[1] - posPV[1]};
    auto momXY = std::array{mom[0], mom[1]};
    auto cos = dotProd(lineDecay, momXY) / std::sqrt(mag2(lineDecay) * mag2(momXY));
    if (cos < -1.) {
      return -1.;
    }
    if (cos > 1.) {
      return 1.;
    }
    return cos;
  }

  /// Calculates proper lifetime times c.
  /// \note Promotes numbers to double before squaring to avoid precision loss in float multiplication.
  /// \param mom  3-momentum array
  /// \param mass  mass
  /// \param length  decay length
  /// \return proper lifetime times c
  template <typename T, typename U, typename V>
  static double ct(const std::array<T, 3>& mom, U length, V mass)
  {
    // c t = l m c^2/(p c)
    return static_cast<double>(length) * static_cast<double>(mass) / p(mom);
  }

  /// Calculates cosine of θ* (theta star).
  /// \note Implemented for 2 prongs only.
  /// \param arrMom  array of two 3-momentum arrays
  /// \param arrMass  array of two masses (in the same order as arrMom)
  /// \param mTot  assumed mass of mother particle
  /// \param iProng  index of the prong
  /// \return cosine of θ* of the i-th prong under the assumption of the invariant mass
  template <typename T, typename U, typename V>
  static double cosThetaStar(const std::array<std::array<T, 3>, 2>& arrMom, const std::array<U, 2>& arrMass, V mTot, int iProng)
  {
    auto pVecTot = pVec(arrMom[0], arrMom[1]);                                                                             // momentum of the mother particle
    auto pTot = p(pVecTot);                                                                                                // magnitude of the momentum of the mother particle
    auto eTot = e(pTot, mTot);                                                                                             // energy of the mother particle
    auto gamma = eTot / mTot;                                                                                              // γ, Lorentz gamma factor of the mother particle
    auto beta = pTot / eTot;                                                                                               // β, velocity of the mother particle
    auto pStar = std::sqrt(sq(sq(mTot) - sq(arrMass[0]) - sq(arrMass[1])) - sq(2 * arrMass[0] * arrMass[1])) / (2 * mTot); // p*, prong momentum in the rest frame of the mother particle
    // p* = √[(M^2 - m1^2 - m2^2)^2 - 4 m1^2 m2^2]/2M
    // Lorentz transformation of the longitudinal momentum of the prong into the detector frame:
    // p_L,i = γ (p*_L,i + β E*_i)
    // p*_L,i = p_L,i/γ - β E*_i
    // cos(θ*_i) = (p_L,i/γ - β E*_i)/p*
    return (dotProd(arrMom[iProng], pVecTot) / (pTot * gamma) - beta * e(pStar, arrMass[iProng])) / pStar;
  }

  /// Sums 3-momenta.
  /// \param args  pack of 3-momentum arrays
  /// \return total 3-momentum array
  template <typename... T>
  static auto pVec(const std::array<T, 3>&... args)
  {
    return sumOfVec(args...);
  }

  /// Calculates momentum squared from momentum components.
  /// \param px,py,pz  {x, y, z} momentum components
  /// \return momentum squared
  static double p2(double px, double py, double pz)
  {
    return sumOfSquares(px, py, pz);
  }

  /// Calculates total momentum squared of a sum of 3-momenta.
  /// \param args  pack of 3-momentum arrays
  /// \return total momentum squared
  template <typename... T>
  static double p2(const std::array<T, 3>&... args)
  {
    return sumOfSquares(getElement(0, args...), getElement(1, args...), getElement(2, args...));
  }

  /// Calculates (total) momentum magnitude.
  /// \param args  {x, y, z} momentum components or pack of 3-momentum arrays
  /// \return (total) momentum magnitude
  template <typename... T>
  static double p(const T&... args)
  {
    return std::sqrt(p2(args...));
  }

  /// Calculates transverse momentum squared from momentum components.
  /// \param px,py  {x, y} momentum components
  /// \return transverse momentum squared
  static double pt2(double px, double py)
  {
    return sumOfSquares(px, py);
  }

  /// Calculates total transverse momentum squared of a sum of 3-(or 2-)momenta.
  /// \param args  pack of 3-(or 2-)momentum arrays
  /// \return total transverse momentum squared
  template <std::size_t N, typename... T>
  static double pt2(const std::array<T, N>&... args)
  {
    return sumOfSquares(getElement(0, args...), getElement(1, args...));
  }

  /// Calculates (total) transverse momentum.
  /// \param args  {x, y} momentum components or pack of 3-(or 2-)momentum arrays
  /// \return (total) transverse momentum
  template <typename... T>
  static double pt(const T&... args)
  {
    return std::sqrt(pt2(args...));
  }

  /// Calculates energy squared from momentum and mass.
  /// \param args  momentum magnitude, mass
  /// \param args  {x, y, z} momentum components, mass
  /// \return energy squared
  template <typename... T>
  static double e2(T... args)
  {
    return sumOfSquares(args...);
  }

  /// Calculates energy squared from momentum vector and mass.
  /// \param mom  3-momentum array
  /// \param mass  mass
  /// \return energy squared
  template <typename T, typename U>
  static double e2(const std::array<T, 3>& mom, U mass)
  {
    return e2(mom[0], mom[1], mom[2], mass);
  }

  /// Calculates energy from momentum and mass.
  /// \param args  momentum magnitude, mass
  /// \param args  {x, y, z} momentum components, mass
  /// \param args  3-momentum array, mass
  /// \return energy
  template <typename... T>
  static double e(const T&... args)
  {
    return std::sqrt(e2(args...));
  }

  /// Calculates invariant mass squared from momentum magnitude and energy.
  /// \param mom  momentum magnitude
  /// \param energy  energy
  /// \return invariant mass squared
  static double m2(double mom, double energy)
  {
    return energy * energy - mom * mom;
  }

  /// Calculates invariant mass squared from momentum aray and energy.
  /// \param mom  3-momentum array
  /// \param energy  energy
  /// \return invariant mass squared
  template <typename T>
  static double m2(const std::array<T, 3>& mom, double energy)
  {
    return energy * energy - p2(mom);
  }

  /// Calculates invariant mass squared from momenta and masses of several particles (prongs).
  /// \param N  number of prongs
  /// \param arrMom  array of N 3-momentum arrays
  /// \param arrMass  array of N masses (in the same order as arrMom)
  /// \return invariant mass squared
  template <std::size_t N, typename T, typename U>
  static double m2(const std::array<std::array<T, 3>, N>& arrMom, const std::array<U, N>& arrMass)
  {
    std::array<double, 3> momTotal{0., 0., 0.}; // candidate momentum vector
    double energyTot{0.};                  // candidate energy
    for (std::size_t iProng = 0; iProng < N; ++iProng) {
      for (std::size_t iMom = 0; iMom < 3; ++iMom) {
        momTotal[iMom] += arrMom[iProng][iMom];
      } // loop over momentum components
      energyTot += e(arrMom[iProng], arrMass[iProng]);
    } // loop over prongs
    return energyTot * energyTot - p2(momTotal);
  }

  /// Calculates invariant mass.
  /// \param args  momentum magnitude, energy
  /// \param args  3-momentum array, energy
  /// \param args  array of momenta, array of masses
  /// \return invariant mass
  template <typename... T>
  static double m(const T&... args)
  {
    return std::sqrt(m2(args...));
  }

  // Calculation of topological quantities

  /// Calculates impact parameter in the bending plane of the particle w.r.t. a point
  /// \param point  {x, y, z} position of the point
  /// \param posSV  {x, y, z} position of the secondary vertex
  /// \param mom  {x, y, z} particle momentum array
  /// \return impact parameter in {x, y}
  template <typename T, typename U, typename V>
  static double impParXY(const T& point, const U& posSV, const std::array<V, 3>& mom)
  {
    // Ported from AliAODRecoDecay::ImpParXY
    auto flightLineXY = std::array{posSV[0] - point[0], posSV[1] - point[1]};
    auto k = dotProd(flightLineXY, std::array{mom[0], mom[1]}) / pt2(mom);
    auto dx = flightLineXY[0] - k * static_cast<double>(mom[0]);
    auto dy = flightLineXY[1] - k * static_cast<double>(mom[1]);
    auto absImpPar = sqrtSumOfSquares(dx, dy);
    auto flightLine = std::array{posSV[0] - point[0], posSV[1] - point[1], posSV[2] - point[2]};
    auto cross = crossProd(mom, flightLine);
    return (cross[2] > 0. ? absImpPar : -1. * absImpPar);
  }

  /// Calculates the difference between measured and expected track impact parameter
  /// normalized to its uncertainty
  /// \param decLenXY decay length in {x, y} plane
  /// \param errDecLenXY error on decay length in {x, y} plane
  /// \param momMother {x, y, z} or {x, y} candidate momentum array
  /// \param impParProng prong impact parameter
  /// \param errImpParProng error on prong impact parameter
  /// \param momProng {x, y, z} or {x, y} prong momentum array
  /// \return normalized difference between expected and observed impact parameter
  template <std::size_t N, std::size_t M, typename T, typename U, typename V, typename W, typename X, typename Y>
  static double normImpParMeasMinusExpProng(T decLenXY, U errDecLenXY, const std::array<V, N>& momMother, W impParProng,
                                            X errImpParProng, const std::array<Y, M>& momProng)
  {
    // Ported from AliAODRecoDecayHF::Getd0MeasMinusExpProng adding normalization directly in the function
    auto sinThetaP = (static_cast<double>(momProng[0]) * static_cast<double>(momMother[1]) - static_cast<double>(momProng[1]) * static_cast<double>(momMother[0])) /
                     (pt(momProng) * pt(momMother));
    auto diff = impParProng - static_cast<double>(decLenXY) * sinThetaP;
    auto errImpParExpProng = static_cast<double>(errDecLenXY) * sinThetaP;
    auto errDiff = sqrtSumOfSquares(errImpParProng, errImpParExpProng);
    return (errDiff > 0. ? diff / errDiff : 0.);
  }

  /// Calculates maximum normalized difference between measured and expected impact parameter of candidate prongs
  /// \param posPV {x, y, z} or {x, y} position of primary vertex
  /// \param posSV {x, y, z} or {x, y} position of secondary vertex
  /// \param errDecLenXY error on decay length in {x, y} plane
  /// \param momMother {x, y, z} or {x, y} candidate momentum array
  /// \param arrImpPar array of prong impact parameters
  /// \param arrErrImpPar array of errors on prong impact parameter (same order as arrImpPar)
  /// \param momMom array of {x, y, z} or {x, y} prong momenta (same order as arrImpPar)
  /// \return maximum normalized difference between expected and observed impact parameters
  template <std::size_t N, std::size_t M, std::size_t K, typename T, typename U, typename V, typename W, typename X,
            typename Y, typename Z>
  static double maxNormalisedDeltaIP(const T& posPV, const U& posSV, V errDecLenXY, const std::array<W, M>& momMother,
                                     const std::array<X, N>& arrImpPar, const std::array<Y, N>& arrErrImpPar,
                                     const std::array<std::array<Z, K>, N>& arrMom)
  {
    auto decLenXY = distanceXY(posPV, posSV);
    double maxNormDeltaIP{0.};
    for (std::size_t iProng = 0; iProng < N; ++iProng) {
      auto prongNormDeltaIP = normImpParMeasMinusExpProng(decLenXY, errDecLenXY, momMother, arrImpPar[iProng],
                                                          arrErrImpPar[iProng], arrMom[iProng]);
      if (std::abs(prongNormDeltaIP) > std::abs(maxNormDeltaIP)) {
        maxNormDeltaIP = prongNormDeltaIP;
      }
    }
    return maxNormDeltaIP;
  }

  /// Finds the mother of an MC particle by looking for the expected PDG code in the mother chain.
  /// \param particlesMC  table with MC particles
  /// \param particle  MC particle
  /// \param PDGMother  expected mother PDG code
  /// \param acceptAntiParticles  switch to accept the antiparticle of the expected mother
  /// \param sign  antiparticle indicator of the found mother w.r.t. PDGMother; 1 if particle, -1 if antiparticle, 0 if mother not found
  /// \param depthMax  maximum decay tree level to check; Mothers up to this level will be considered. If -1, all levels are considered.
  /// \return index of the mother particle if found, -1 otherwise
  template <bool acceptFlavourOscillation = false, typename T>
  static int getMother(const T& particlesMC,
                       const typename T::iterator& particle,
                       int PDGMother,
                       bool acceptAntiParticles = false,
                       int8_t* sign = nullptr,
                       int8_t depthMax = -1)
  {
    int8_t sgn = 0;           // 1 if the expected mother is particle, -1 if antiparticle (w.r.t. PDGMother)
    int indexMother = -1;     // index of the final matched mother, if found
    int stage = 0;            // mother tree level (just for debugging)
    bool motherFound = false; // true when the desired mother particle is found in the kine tree
    if (sign) {
      *sign = sgn;
    }

    // vector of vectors with mother indices; each line corresponds to a "stage"
    std::vector<std::vector<int64_t>> arrayIds{};
    std::vector<int64_t> initVec{particle.globalIndex()};
    arrayIds.push_back(initVec); // the first vector contains the index of the original particle

    while (!motherFound && arrayIds[-stage].size() > 0 && (depthMax < 0 || -stage < depthMax)) {
      // vector of mother indices for the current stage
      std::vector<int64_t> arrayIdsStage{};
      for (auto& iPart : arrayIds[-stage]) { // check all the particles that were the mothers at the previous stage
        auto particleMother = particlesMC.rawIteratorAt(iPart - particlesMC.offset());
        if (particleMother.has_mothers()) {
          for (auto iMother = particleMother.mothersIds().front(); iMother <= particleMother.mothersIds().back(); ++iMother) { // loop over the mother particles of the analysed particle
            if (std::find(arrayIdsStage.begin(), arrayIdsStage.end(), iMother) != arrayIdsStage.end()) {                       // if a mother is still present in the vector, do not check it again
              continue;
            }
            auto mother = particlesMC.rawIteratorAt(iMother - particlesMC.offset());
            // Check mother's PDG code.
            auto PDGParticleIMother = mother.pdgCode(); // PDG code of the mother
            // printf("getMother: ");
            // for (int i = stage; i < 0; i++) // Indent to make the tree look nice.
            //   printf(" ");
            // printf("Stage %d: Mother PDG: %d, Index: %d\n", stage, PDGParticleIMother, iMother);
            if (PDGParticleIMother == PDGMother) { // exact PDG match
              sgn = 1;
              indexMother = iMother;
              motherFound = true;
              break;
            } else if (acceptAntiParticles && PDGParticleIMother == -PDGMother) { // antiparticle PDG match
              sgn = -1;
              indexMother = iMother;
              motherFound = true;
              break;
            }
            // add mother index in the vector for the current stage
            arrayIdsStage.push_back(iMother);
          }
        }
      }
      // add vector of mother indices for the current stage
      arrayIds.push_back(arrayIdsStage);
      stage--;
    }
    if (sign) {
      if constexpr (acceptFlavourOscillation) {
        if (std::abs(particle.getGenStatusCode()) == PdgStatusCodeAfterFlavourOscillation) { // take possible flavour oscillation of B0(s) mother into account
          sgn *= -1;                                                                         // select the sign of the mother after oscillation (and not before)
        }
      }
      *sign = sgn;
    }

    return indexMother;
  }

  /// Gets the complete list of indices of final-state daughters of an MC particle.
  /// \param particle  MC particle
  /// \param list  vector where the indices of final-state daughters will be added
  /// \param arrPDGFinal  array of PDG codes of particles to be considered final if found
  /// \param depthMax  maximum decay tree level; Daughters at this level (or beyond) will be considered final. If -1, all levels are considered.
  /// \param stage  decay tree level; If different from 0, the particle itself will be added in the list in case it has no daughters.
  /// \note Final state is defined as particles from arrPDGFinal plus final daughters of any other decay branch.
  /// \note Antiparticles of particles in arrPDGFinal are accepted as well.
  template <std::size_t N, typename T>
  static void getDaughters(const T& particle,
                           std::vector<int>* list,
                           const std::array<int, N>& arrPDGFinal,
                           int8_t depthMax = -1,
                           int8_t stage = 0)
  {
    if (!list) {
      // Printf("getDaughters: Error: No list!");
      return;
    }
    bool isFinal = false;                     // Flag to indicate the end of recursion
    if (depthMax > -1 && stage >= depthMax) { // Maximum depth has been reached (or exceeded).
      isFinal = true;
    }
    // Check whether there are any daughters.
    if (!isFinal && !particle.has_daughters()) {
      // If the original particle has no daughters, we do nothing and exit.
      if (stage == 0) {
        // Printf("getDaughters: No daughters of %d", index);
        return;
      }
      // If this is not the original particle, we are at the end of this branch and this particle is final.
      isFinal = true;
    }
    auto PDGParticle = std::abs(particle.pdgCode());
    // If this is not the original particle, check its PDG code.
    if (!isFinal && stage > 0) {
      // If the particle has daughters but is considered to be final, we label it as final.
      for (auto PDGi : arrPDGFinal) {
        if (PDGParticle == std::abs(PDGi)) { // Accept antiparticles.
          isFinal = true;
          break;
        }
      }
    }
    // If the particle is labelled as final, we add this particle in the list of final daughters and exit.
    if (isFinal) {
      // printf("getDaughters: ");
      // for (int i = 0; i < stage; i++) // Indent to make the tree look nice.
      //   printf(" ");
      // printf("Stage %d: Adding %d (PDG %d) as final daughter.\n", stage, index, PDGParticle);
      list->push_back(particle.globalIndex());
      return;
    }
    // If we are here, we have to follow the daughter tree.
    // printf("getDaughters: ");
    // for (int i = 0; i < stage; i++) // Indent to make the tree look nice.
    //  printf(" ");
    // printf("Stage %d: %d (PDG %d) -> %d-%d\n", stage, index, PDGParticle, indexDaughterFirst, indexDaughterLast);
    // Call itself to get daughters of daughters recursively.
    stage++;
    for (auto& dau : particle.template daughters_as<typename std::decay_t<T>::parent_t>()) {
      getDaughters(dau, list, arrPDGFinal, depthMax, stage);
    }
  }

  /// Checks whether the reconstructed decay candidate is the expected decay.
  /// \param particlesMC  table with MC particles
  /// \param arrDaughters  array of candidate daughters
  /// \param PDGMother  expected mother PDG code
  /// \param arrPDGDaughters  array of expected daughter PDG codes
  /// \param acceptAntiParticles  switch to accept the antiparticle version of the expected decay
  /// \param sign  antiparticle indicator of the found mother w.r.t. PDGMother; 1 if particle, -1 if antiparticle, 0 if mother not found
  /// \param depthMax  maximum decay tree level to check; Daughters up to this level will be considered. If -1, all levels are considered.
  /// \return index of the mother particle if the mother and daughters are correct, -1 otherwise
  template <bool acceptFlavourOscillation = false, std::size_t N, typename T, typename U>
  static int getMatchedMCRec(const T& particlesMC,
                             const std::array<U, N>& arrDaughters,
                             int PDGMother,
                             std::array<int, N> arrPDGDaughters,
                             bool acceptAntiParticles = false,
                             int8_t* sign = nullptr,
                             int depthMax = 1)
  {
    // Printf("MC Rec: Expected mother PDG: %d", PDGMother);
    int8_t coefFlavourOscillation = 1;     // 1 if no B0(s) flavour oscillation occured, -1 else
    int8_t sgn = 0;                        // 1 if the expected mother is particle, -1 if antiparticle (w.r.t. PDGMother)
    int indexMother = -1;                  // index of the mother particle
    std::vector<int> arrAllDaughtersIndex; // vector of indices of all daughters of the mother of the first provided daughter
    std::array<int, N> arrDaughtersIndex;  // array of indices of provided daughters
    if (sign) {
      *sign = sgn;
    }
    if constexpr (acceptFlavourOscillation) {
      // Loop over decay candidate prongs to spot possible oscillation decay product
      for (std::size_t iProng = 0; iProng < N; ++iProng) {
        if (!arrDaughters[iProng].has_mcParticle()) {
          return -1;
        }
        auto particleI = arrDaughters[iProng].mcParticle();                                   // ith daughter particle
        if (std::abs(particleI.getGenStatusCode()) == PdgStatusCodeAfterFlavourOscillation) { // oscillation decay product spotted
          coefFlavourOscillation = -1;                                                        // select the sign of the mother after oscillation (and not before)
          break;
        }
      }
    }
    // Loop over decay candidate prongs
    for (std::size_t iProng = 0; iProng < N; ++iProng) {
      if (!arrDaughters[iProng].has_mcParticle()) {
        return -1;
      }
      auto particleI = arrDaughters[iProng].mcParticle(); // ith daughter particle
      arrDaughtersIndex[iProng] = particleI.globalIndex();
      // Get the list of daughter indices from the mother of the first prong.
      if (iProng == 0) {
        // Get the mother index and its sign.
        // PDG code of the first daughter's mother determines whether the expected mother is a particle or antiparticle.
        indexMother = getMother(particlesMC, particleI, PDGMother, acceptAntiParticles, &sgn, depthMax);
        // Check whether mother was found.
        if (indexMother <= -1) {
          // Printf("MC Rec: Rejected: bad mother index or PDG");
          return -1;
        }
        // Printf("MC Rec: Good mother: %d", indexMother);
        auto particleMother = particlesMC.rawIteratorAt(indexMother - particlesMC.offset());
        // Check the daughter indices.
        if (!particleMother.has_daughters()) {
          // Printf("MC Rec: Rejected: bad daughter index range: %d-%d", particleMother.daughtersIds().front(), particleMother.daughtersIds().back());
          return -1;
        }
        // Check that the number of direct daughters is not larger than the number of expected final daughters.
        if (particleMother.daughtersIds().back() - particleMother.daughtersIds().front() + 1 > static_cast<int>(N)) {
          // Printf("MC Rec: Rejected: too many direct daughters: %d (expected %ld final)", particleMother.daughtersIds().back() - particleMother.daughtersIds().front() + 1, N);
          return -1;
        }
        // Get the list of actual final daughters.
        getDaughters(particleMother, &arrAllDaughtersIndex, arrPDGDaughters, depthMax);
        // printf("MC Rec: Mother %d has %d final daughters:", indexMother, arrAllDaughtersIndex.size());
        // for (auto i : arrAllDaughtersIndex) {
        //   printf(" %d", i);
        // }
        // printf("\n");
        //  Check whether the number of actual final daughters is equal to the number of expected final daughters (i.e. the number of provided prongs).
        if (arrAllDaughtersIndex.size() != N) {
          // Printf("MC Rec: Rejected: incorrect number of final daughters: %ld (expected %ld)", arrAllDaughtersIndex.size(), N);
          return -1;
        }
      }
      // Check that the daughter is in the list of final daughters.
      // (Check that the daughter is not a stepdaughter, i.e. particle pointing to the mother while not being its daughter.)
      bool isDaughterFound = false; // Is the index of this prong among the remaining expected indices of daughters?
      for (std::size_t iD = 0; iD < arrAllDaughtersIndex.size(); ++iD) {
        if (arrDaughtersIndex[iProng] == arrAllDaughtersIndex[iD]) {
          arrAllDaughtersIndex[iD] = -1; // Remove this index from the array of expected daughters. (Rejects twin daughters, i.e. particle considered twice as a daughter.)
          isDaughterFound = true;
          break;
        }
      }
      if (!isDaughterFound) {
        // Printf("MC Rec: Rejected: bad daughter index: %d not in the list of final daughters", arrDaughtersIndex[iProng]);
        return -1;
      }
      // Check daughter's PDG code.
      auto PDGParticleI = particleI.pdgCode(); // PDG code of the ith daughter
      // Printf("MC Rec: Daughter %d PDG: %d", iProng, PDGParticleI);
      bool isPDGFound = false; // Is the PDG code of this daughter among the remaining expected PDG codes?
      for (std::size_t iProngCp = 0; iProngCp < N; ++iProngCp) {
        if (PDGParticleI == coefFlavourOscillation * sgn * arrPDGDaughters[iProngCp]) {
          arrPDGDaughters[iProngCp] = 0; // Remove this PDG code from the array of expected ones.
          isPDGFound = true;
          break;
        }
      }
      if (!isPDGFound) {
        // Printf("MC Rec: Rejected: bad daughter PDG: %d", PDGParticleI);
        return -1;
      }
    }
    // Printf("MC Rec: Accepted: m: %d", indexMother);
    if (sign) {
      *sign = sgn;
    }
    return indexMother;
  }

  /// Checks whether the MC particle is the expected one.
  /// \param particlesMC  table with MC particles
  /// \param candidate  candidate MC particle
  /// \param PDGParticle  expected particle PDG code
  /// \param acceptAntiParticles  switch to accept the antiparticle
  /// \param sign  antiparticle indicator of the candidate w.r.t. PDGParticle; 1 if particle, -1 if antiparticle, 0 if not matched
  /// \return true if PDG code of the particle is correct, false otherwise
  template <bool acceptFlavourOscillation = false, typename T, typename U>
  static int isMatchedMCGen(const T& particlesMC,
                            const U& candidate,
                            int PDGParticle,
                            bool acceptAntiParticles = false,
                            int8_t* sign = nullptr)
  {
    std::array<int, 0> arrPDGDaughters;
    return isMatchedMCGen<acceptFlavourOscillation>(particlesMC, candidate, PDGParticle, std::move(arrPDGDaughters), acceptAntiParticles, sign);
  }

  /// Check whether the MC particle is the expected one and whether it decayed via the expected decay channel.
  /// \param particlesMC  table with MC particles
  /// \param candidate  candidate MC particle
  /// \param PDGParticle  expected particle PDG code
  /// \param arrPDGDaughters  array of expected PDG codes of daughters
  /// \param acceptAntiParticles  switch to accept the antiparticle
  /// \param sign  antiparticle indicator of the candidate w.r.t. PDGParticle; 1 if particle, -1 if antiparticle, 0 if not matched
  /// \param depthMax  maximum decay tree level to check; Daughters up to this level will be considered. If -1, all levels are considered.
  /// \param listIndexDaughters  vector of indices of found daughter
  /// \return true if PDG codes of the particle and its daughters are correct, false otherwise
  template <bool acceptFlavourOscillation = false, std::size_t N, typename T, typename U>
  static bool isMatchedMCGen(const T& particlesMC,
                             const U& candidate,
                             int PDGParticle,
                             std::array<int, N> arrPDGDaughters,
                             bool acceptAntiParticles = false,
                             int8_t* sign = nullptr,
                             int depthMax = 1,
                             std::vector<int>* listIndexDaughters = nullptr)
  {
    // Printf("MC Gen: Expected particle PDG: %d", PDGParticle);
    int8_t coefFlavourOscillation = 1; // 1 if no B0(s) flavour oscillation occured, -1 else
    int8_t sgn = 0; // 1 if the expected mother is particle, -1 if antiparticle (w.r.t. PDGParticle)
    if (sign) {
      *sign = sgn;
    }
    // Check the PDG code of the particle.
    auto PDGCandidate = candidate.pdgCode();
    // Printf("MC Gen: Candidate PDG: %d", PDGCandidate);
    if (PDGCandidate == PDGParticle) { // exact PDG match
      sgn = 1;
    } else if (acceptAntiParticles && PDGCandidate == -PDGParticle) { // antiparticle PDG match
      sgn = -1;
    } else {
      // Printf("MC Gen: Rejected: bad particle PDG: %s%d != %d", acceptAntiParticles ? "abs " : "", PDGCandidate, std::abs(PDGParticle));
      return false;
    }
    // Check the PDG codes of the decay products.
    if (N > 0) {
      // Printf("MC Gen: Checking %d daughters", N);
      std::vector<int> arrAllDaughtersIndex; // vector of indices of all daughters
      // Check the daughter indices.
      if (!candidate.has_daughters()) {
        // Printf("MC Gen: Rejected: bad daughter index range: %d-%d", candidate.daughtersIds().front(), candidate.daughtersIds().back());
        return false;
      }
      // Check that the number of direct daughters is not larger than the number of expected final daughters.
      if (candidate.daughtersIds().back() - candidate.daughtersIds().front() + 1 > static_cast<int>(N)) {
        // Printf("MC Gen: Rejected: too many direct daughters: %d (expected %ld final)", candidate.daughtersIds().back() - candidate.daughtersIds().front() + 1, N);
        return false;
      }
      // Get the list of actual final daughters.
      getDaughters(candidate, &arrAllDaughtersIndex, arrPDGDaughters, depthMax);
      // printf("MC Gen: Mother %ld has %ld final daughters:", candidate.globalIndex(), arrAllDaughtersIndex.size());
      // for (auto i : arrAllDaughtersIndex) {
      //   printf(" %d", i);
      // }
      // printf("\n");
      //  Check whether the number of final daughters is equal to the required number.
      if (arrAllDaughtersIndex.size() != N) {
        // Printf("MC Gen: Rejected: incorrect number of final daughters: %ld (expected %ld)", arrAllDaughtersIndex.size(), N);
        return false;
      }
      if constexpr (acceptFlavourOscillation) {
        // Loop over decay candidate prongs to spot possible oscillation decay product
        for (auto indexDaughterI : arrAllDaughtersIndex) {
          auto candidateDaughterI = particlesMC.rawIteratorAt(indexDaughterI - particlesMC.offset());    // ith daughter particle
          if (std::abs(candidateDaughterI.getGenStatusCode()) == PdgStatusCodeAfterFlavourOscillation) { // oscillation decay product spotted
            coefFlavourOscillation = -1;                                                                 // select the sign of the mother after oscillation (and not before)
            break;
          }
        }
      }
      // Check daughters' PDG codes.
      for (auto indexDaughterI : arrAllDaughtersIndex) {
        auto candidateDaughterI = particlesMC.rawIteratorAt(indexDaughterI - particlesMC.offset()); // ith daughter particle
        auto PDGCandidateDaughterI = candidateDaughterI.pdgCode();                                  // PDG code of the ith daughter
        // Printf("MC Gen: Daughter %d PDG: %d", indexDaughterI, PDGCandidateDaughterI);
        bool isPDGFound = false; // Is the PDG code of this daughter among the remaining expected PDG codes?
        for (std::size_t iProngCp = 0; iProngCp < N; ++iProngCp) {
          if (PDGCandidateDaughterI == coefFlavourOscillation * sgn * arrPDGDaughters[iProngCp]) {
            arrPDGDaughters[iProngCp] = 0; // Remove this PDG code from the array of expected ones.
            isPDGFound = true;
            break;
          }
        }
        if (!isPDGFound) {
          // Printf("MC Gen: Rejected: bad daughter PDG: %d", PDGCandidateDaughterI);
          return false;
        }
      }
      if (listIndexDaughters) {
        *listIndexDaughters = arrAllDaughtersIndex;
      }
    }
    // Printf("MC Gen: Accepted: m: %d", candidate.globalIndex());
    if (sign) {
      *sign = sgn;
    }
    return true;
  }

  /// Finds the origin (from charm hadronisation or beauty-hadron decay) of charm hadrons. It can be used also to verify whether a particle derives from a charm or beauty decay.
  /// \param particlesMC  table with MC particles
  /// \param particle  MC particle
  /// \param searchUpToQuark if true tag origin based on charm/beauty quark otherwise on the presence of a b-hadron or c-hadron, with c-hadrons themselves marked as prompt
  /// \return an integer corresponding to the origin (0: none, 1: prompt, 2: nonprompt) as in OriginType
  template <typename T>
  static int getCharmHadronOrigin(const T& particlesMC,
                                  const typename T::iterator& particle,
                                  const bool searchUpToQuark = false)
  {
    int stage = 0; // mother tree level (just for debugging)

    // vector of vectors with mother indices; each line corresponds to a "stage"
    std::vector<std::vector<int64_t>> arrayIds{};
    std::vector<int64_t> initVec{particle.globalIndex()};
    arrayIds.push_back(initVec); // the first vector contains the index of the original particle
    auto PDGParticle = std::abs(particle.pdgCode());
    bool couldBePrompt = false;
    if (PDGParticle / 100 == 4 || PDGParticle / 1000 == 4) {
      couldBePrompt = true;
    }
    while (arrayIds[-stage].size() > 0) {
      // vector of mother indices for the current stage
      std::vector<int64_t> arrayIdsStage{};
      for (auto& iPart : arrayIds[-stage]) { // check all the particles that were the mothers at the previous stage
        auto particleMother = particlesMC.rawIteratorAt(iPart - particlesMC.offset());
        if (particleMother.has_mothers()) {
          for (auto iMother = particleMother.mothersIds().front(); iMother <= particleMother.mothersIds().back(); ++iMother) { // loop over the mother particles of the analysed particle
            if (std::find(arrayIdsStage.begin(), arrayIdsStage.end(), iMother) != arrayIdsStage.end()) {                       // if a mother is still present in the vector, do not check it again
              continue;
            }
            auto mother = particlesMC.rawIteratorAt(iMother - particlesMC.offset());
            // Check mother's PDG code.
            auto PDGParticleIMother = std::abs(mother.pdgCode()); // PDG code of the mother
            // printf("getMother: ");
            // for (int i = stage; i < 0; i++) // Indent to make the tree look nice.
            //   printf(" ");
            // printf("Stage %d: Mother PDG: %d, Index: %d\n", stage, PDGParticleIMother, iMother);

            if (searchUpToQuark) {
              if (PDGParticleIMother == 5) { // b quark
                return OriginType::NonPrompt;
              }
              if (PDGParticleIMother == 4) { // c quark
                return OriginType::Prompt;
              }
            } else {
              if (
                (PDGParticleIMother / 100 == 5 || // b mesons
                 PDGParticleIMother / 1000 == 5)  // b baryons
              ) {
                return OriginType::NonPrompt;
              }
              if (
                (PDGParticleIMother / 100 == 4 || // c mesons
                 PDGParticleIMother / 1000 == 4)  // c baryons
              ) {
                couldBePrompt = true;
              }
            }
            // add mother index in the vector for the current stage
            arrayIdsStage.push_back(iMother);
          }
        }
      }
      // add vector of mother indices for the current stage
      arrayIds.push_back(arrayIdsStage);
      stage--;
    }
    if (!searchUpToQuark && couldBePrompt) { // Returns prompt if it's a charm hadron or a charm-hadron daughter. Note: 1) LF decay particles from cases like -> Lc -> p K0S, K0S -> pi pi are marked as prompt. 2) if particles from HF parton showers have to be searched, switch to option "search up to quark"
      return OriginType::Prompt;
    }
    return OriginType::None;
  }
};

#endif // COMMON_CORE_RECODECAY_H_
