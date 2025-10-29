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

/// \file FemtoUniverseSpherHarMath.h
/// \brief Definition of the FemtoUniverseMath Container for the calculations of spherical harmonics components
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESPHERHARMATH_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESPHERHARMATH_H_

#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <algorithm>
#include <complex>
#include <vector>

namespace o2::analysis::femto_universe
{

/// \class FemtoUniverseMath
/// \brief Container for math calculations of quantities related to pairs
class FemtoUniverseSpherHarMath
{
 public:
  /// Values of various coefficients
  void initializeYlms()
  {
    double oneoversqrtpi = 1.0 / std::sqrt(o2::constants::math::PI);

    // l=0 prefactors
    fgPrefactors[0] = 0.5 * oneoversqrtpi;

    // l=1 prefactors
    fgPrefactors[1] = 0.5 * std::sqrt(3.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[2] = 0.5 * std::sqrt(3.0) * oneoversqrtpi;
    fgPrefactors[3] = -fgPrefactors[1];

    // l=2 prefactors
    fgPrefactors[4] = 0.25 * std::sqrt(15.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[5] = 0.5 * std::sqrt(15.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[6] = 0.25 * std::sqrt(5.0) * oneoversqrtpi;
    fgPrefactors[7] = -fgPrefactors[5];
    fgPrefactors[8] = fgPrefactors[4];

    // l=3 prefactors
    fgPrefactors[9] = 0.125 * std::sqrt(35.0) * oneoversqrtpi;
    fgPrefactors[10] = 0.25 * std::sqrt(105.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[11] = 0.125 * std::sqrt(21.0) * oneoversqrtpi;
    fgPrefactors[12] = 0.25 * std::sqrt(7.0) * oneoversqrtpi;
    fgPrefactors[13] = -fgPrefactors[11];
    fgPrefactors[14] = fgPrefactors[10];
    fgPrefactors[15] = -fgPrefactors[9];

    // l=4 prefactors
    fgPrefactors[16] = 3.0 / 16.0 * std::sqrt(35.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[17] = 3.0 / 8.0 * std::sqrt(35.0) * oneoversqrtpi;
    fgPrefactors[18] = 3.0 / 8.0 * std::sqrt(5.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[19] = 3.0 / 8.0 * std::sqrt(5.0) * oneoversqrtpi;
    fgPrefactors[20] = 3.0 / 16.0 * oneoversqrtpi;
    fgPrefactors[21] = -fgPrefactors[19];
    fgPrefactors[22] = fgPrefactors[18];
    fgPrefactors[23] = -fgPrefactors[17];
    fgPrefactors[24] = fgPrefactors[16];

    // l=5 prefactors
    fgPrefactors[25] = 3.0 / 32.0 * std::sqrt(77.0) * oneoversqrtpi;
    fgPrefactors[26] = 3.0 / 16.0 * std::sqrt(385.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[27] = 1.0 / 32.0 * std::sqrt(385.0) * oneoversqrtpi;
    fgPrefactors[28] = 1.0 / 8.0 * std::sqrt(1155.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[29] = 1.0 / 16.0 * std::sqrt(165.0 / 2.0) * oneoversqrtpi;
    fgPrefactors[30] = 1.0 / 16.0 * std::sqrt(11.0) * oneoversqrtpi;
    fgPrefactors[31] = -fgPrefactors[29];
    fgPrefactors[32] = fgPrefactors[28];
    fgPrefactors[33] = -fgPrefactors[27];
    fgPrefactors[34] = fgPrefactors[26];
    fgPrefactors[35] = -fgPrefactors[25];

    fgPrefshift[0] = 0;
    fgPrefshift[1] = 2;
    fgPrefshift[2] = 6;
    fgPrefshift[3] = 12;
    fgPrefshift[4] = 20;
    fgPrefshift[5] = 30;

    fgPlmshift[0] = 0;
    fgPlmshift[1] = 2;
    fgPlmshift[2] = 5;
    fgPlmshift[3] = 9;
    fgPlmshift[4] = 14;
    fgPlmshift[5] = 20;
  }

  /// Function to calculate Legendre Polynomials
  ///  \param lmax Maximum value of L component
  ///  \param ctheta Value of theta
  ///  \param lbuf values of coefficients
  void legendreUpToYlm(int lmax, double ctheta, double* lbuf)
  {
    // Calculate a set of legendre polynomials up to a given l
    // with spherical input
    double sins[6];
    double coss[6];
    sins[0] = 0.0;
    coss[0] = 1.0;
    sins[1] = std::sqrt(1 - ctheta * ctheta);
    coss[1] = ctheta;
    for (int iter = 2; iter < 6; iter++) {
      sins[iter] = sins[iter - 1] * sins[1];
      coss[iter] = coss[iter - 1] * coss[1];
    }

    // Legendre polynomials l=0
    lbuf[0] = 1.0;

    // Legendre polynomials l=1
    if (lmax > 0) {
      lbuf[1] = sins[1];
      lbuf[2] = coss[1];
    }

    // Legendre polynomials l=2
    if (lmax > 1) {
      lbuf[3] = sins[2];
      lbuf[4] = sins[1] * coss[1];
      lbuf[5] = 3 * coss[2] - 1;
    }

    // Legendre polynomials l=3
    if (lmax > 2) {
      lbuf[6] = sins[3];
      lbuf[7] = sins[2] * coss[1];
      lbuf[8] = (5 * coss[2] - 1) * sins[1];
      lbuf[9] = 5 * coss[3] - 3 * coss[1];
    }

    // Legendre polynomials l=4
    if (lmax > 3) {
      lbuf[10] = sins[4];
      lbuf[11] = sins[3] * coss[1];
      lbuf[12] = (7 * coss[2] - 1) * sins[2];
      lbuf[13] = (7 * coss[3] - 3 * coss[1]) * sins[1];
      lbuf[14] = 35 * coss[4] - 30 * coss[2] + 3;
    }

    // Legendre polynomials l=5
    if (lmax > 4) {
      lbuf[15] = sins[5];
      lbuf[16] = sins[4] * coss[1];
      lbuf[17] = (9 * coss[2] - 1) * sins[3];
      lbuf[18] = (3 * coss[3] - 1 * coss[1]) * sins[2];
      lbuf[19] = (21 * coss[4] - 14 * coss[2] + 1) * sins[1];
      lbuf[20] = 63 * coss[5] - 70 * coss[3] + 15 * coss[1];
    }
  }

  /// Function to calculate a set of Ylms up to a given l with cartesian input
  void doYlmUpToL(int lmax, double x, double y, double z, std::complex<double>* ylms)
  {
    double ctheta, phi;

    double r = std::sqrt(x * x + y * y + z * z);
    if (r < 1e-10 || std::fabs(z) < 1e-10)
      ctheta = 0.0;
    else
      ctheta = z / r;
    phi = std::atan2(y, x);
    doYlmUpToL(lmax, ctheta, phi, ylms);
  }

  /// Function to calculate a set of Ylms up to a given l with spherical input
  void doYlmUpToL(int lmax, double ctheta, double phi, std::complex<double>* ylms)
  {
    int lcur = 0;
    double lpol;

    double coss[6];
    double sins[6];

    double lbuf[36];
    legendreUpToYlm(lmax, ctheta, lbuf);
    initializeYlms();

    for (int iter = 1; iter <= lmax; iter++) {
      coss[iter - 1] = std::cos(iter * phi);
      sins[iter - 1] = std::sin(iter * phi);
    }

    ylms[lcur++] = fgPrefactors[0] * lbuf[0] * std::complex<double>(1, 0);

    for (int il = 1; il <= lmax; il++) {
      // First im = 0
      ylms[lcur + il] = fgPrefactors[fgPrefshift[il]] * lbuf[static_cast<int>(fgPlmshift[il])] * std::complex<double>(1.0, 0.0);
      // Im != 0
      for (int im = 1; im <= il; im++) {
        lpol = lbuf[static_cast<int>(fgPlmshift[il]) - im];
        ylms[lcur + il - im] = fgPrefactors[fgPrefshift[il] - im] * lpol * std::complex<double>(coss[im - 1], -sins[im - 1]);
        ylms[lcur + il + im] = fgPrefactors[fgPrefshift[il] + im] * lpol * std::complex<double>(coss[im - 1], sins[im - 1]);
      }
      lcur += 2 * il + 1;
    }
  }

 private:
  static std::complex<double> fCeiphi(double phi);

  std::array<float, 36> fgPrefactors;
  std::array<float, 10> fgPrefshift;
  std::array<float, 10> fgPlmshift;
};

} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESPHERHARMATH_H_
