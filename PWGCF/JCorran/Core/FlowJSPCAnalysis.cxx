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

// Header files.

// O2 headers.

// O2 Physics headers.

#include "PWGCF/JCorran/Core/FlowJSPCAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace std;

TComplex FlowJSPCAnalysis::q(const int harmN, const int p)
{
  if (harmN >= 0)
    return qvecs->QvectorQC[harmN][p];
  return TComplex::Conjugate(qvecs->QvectorQC[-harmN][p]);
} // End of Q

/// \brief Calculate multi-particle correlators using recursion.
/// \param n Number of particles in the correlator.
/// \param harmonic Array of length n of harmonics.
/// \return Complex value of the multiparticle correlator.
/// \note Improved faster version) originally developed by Kristjan Gulbrandsen
/// (gulbrand@nbi.dk).
TComplex FlowJSPCAnalysis::recursion(int n, int* harmonic, int mult = 1, int skip = 0)
{
  int nm1 = n - 1;
  TComplex c(q(harmonic[nm1], mult));
  if (nm1 == 0)
    return c;
  c *= recursion(nm1, harmonic);
  if (nm1 == skip)
    return c;

  int multp1 = mult + 1;
  int nm2 = n - 2;
  int counter1 = 0;
  int hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(recursion(nm1, harmonic, multp1, nm2));
  int counter2 = n - 3;

  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1)
    return c - c2;
  return c - static_cast<double>(mult) * c2;
} // End of recursion

void FlowJSPCAnalysis::calculateCorrelators(const int fCentBin)
{
  // Loop over the combinations of harmonics and calculate the corresponding SPC num and den.

  // Declare the arrays to later fill all the needed bins for the correlators
  // and the error terms.
  double* dataCorrelation = new double[3]; // cosine, weight, sine.
  double correlationNum;
  double weightCorrelationNum;
  double correlationDenom;
  double weightCorrelationDenom;

  for (int i = 0; i < 14; ++i)
    fCorrelDenoms[i] = 0;

  for (int j = 0; j < 12; j++) {
    if (fHarmosArray[j][0] == 0) {
      continue;
    } // Skip null correlator list.

    // Calculate the numerator.
    int hArrayNum[7] = {0};
    for (int iH = 0; iH < 7; iH++) {
      hArrayNum[iH] = fHarmosArray[j][iH + 1];
    }
    correlation(fHarmosArray[j][0], 7, hArrayNum, dataCorrelation);
    correlationNum = dataCorrelation[0];
    weightCorrelationNum = dataCorrelation[1];

    // Calculate the denominator.
    int nPartDen = 2 * fHarmosArray[j][0];
    int hArrayDen[14] = {0};
    for (int iH = 0; iH < 7; iH++) {
      hArrayDen[2 * iH] = hArrayNum[iH];
      hArrayDen[2 * iH + 1] = -1 * hArrayNum[iH];
    }

    correlation(nPartDen, 14, hArrayDen, dataCorrelation);
    correlationDenom = dataCorrelation[0];
    weightCorrelationDenom = dataCorrelation[1];

    // Check if the values are real numbers before filling.
    if (std::isnan(correlationNum) || std::isnan(correlationDenom) || std::isnan(weightCorrelationNum) || std::isnan(weightCorrelationDenom))
      continue;

    // Histogram filling
    fillHistograms(fCentBin, j, correlationNum, correlationDenom, weightCorrelationNum, weightCorrelationDenom);

    correlationNum = 0.;
    weightCorrelationNum = 0.;
    correlationDenom = 0.;
    weightCorrelationDenom = 0.;
  }
}

void FlowJSPCAnalysis::fillHistograms(const int fCentBin, int ind, double cNum, double cDenom, double wNum, double wDenom)
{
  switch (fCentBin) {
    case 0: {
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 1: {
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 2: {
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 3: {
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 4: {
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 5: {
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 6: {
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 7: {
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    case 8: {
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("fResults"), 2. * static_cast<float>(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("fResults"), 2. * static_cast<float>(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("fCovResults"), 4. * static_cast<float>(ind) + 3.5, wDenom, 1.);
    } break;
    default:
      return;
  }
}

void FlowJSPCAnalysis::fillQAHistograms(const int fCentBin, double phi, double phiWeight)
{
  switch (fCentBin) {
    case 0: {
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[0]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 1: {
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[1]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 2: {
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[2]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 3: {
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[3]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 4: {
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[4]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 5: {
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[5]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 6: {
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[6]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 7: {
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[7]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 8: {
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(MCentClasses[8]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    default:
      return;
  }
}

void FlowJSPCAnalysis::correlation(int c_nPart, int c_nHarmo, int* harmo, double* correlData)
{
  // Calculate the correlators for the provided set of harmonics using Q-vectors.
  // Protection against anisotropic correlators.
  int sumHarmo = 0;
  for (int i = 0; i < c_nHarmo; i++) {
    sumHarmo += harmo[i];
  }
  if (sumHarmo != 0) {
    LOGF(error, "\nOups, this correlator is not isotropic(sum = %d). Bye\n", sumHarmo);
    return;
  }

  switch (c_nPart) {
    case 2: {
      int harmonicsTwoNum[2] = {harmo[0], harmo[1]};
      int harmonicsTwoDen[2] = {0, 0};

      if (!fCorrelDenoms[1]) {
        fCorrelDenoms[1] = recursion(2, harmonicsTwoDen).Re();
      }

      TComplex twoRecursion = recursion(2, harmonicsTwoNum) / fCorrelDenoms[1];

      correlData[0] = twoRecursion.Re(); // <cos(h1*phi1+h2*phi2)>
      correlData[1] = fCorrelDenoms[1];  // weight
      correlData[2] = twoRecursion.Im(); // <sin(h1*phi1+h2*phi2)>
    } break;
    case 3: {
      int harmonicsThreeNum[3] = {harmo[0], harmo[1], harmo[2]};
      int harmonicsThreeDen[3] = {0, 0, 0};

      if (!fCorrelDenoms[2]) {
        fCorrelDenoms[2] = recursion(3, harmonicsThreeDen).Re();
      }

      TComplex threeRecursion = recursion(3, harmonicsThreeNum) / fCorrelDenoms[2];

      correlData[0] = threeRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3)>
      correlData[1] = fCorrelDenoms[2];    // weight
      correlData[2] = threeRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3)>
    } break;
    case 4: {
      int harmonicsFourNum[4] = {harmo[0], harmo[1], harmo[2], harmo[3]};
      int harmonicsFourDen[4] = {0, 0, 0, 0};

      if (!fCorrelDenoms[3]) {
        fCorrelDenoms[3] = recursion(4, harmonicsFourDen).Re();
      }

      TComplex fourRecursion = recursion(4, harmonicsFourNum) / fCorrelDenoms[3];

      correlData[0] = fourRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>
      correlData[1] = fCorrelDenoms[3];   // weight
      correlData[2] = fourRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>
    } break;
    case 5: {
      int harmonicsFiveNum[5] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4]};
      int harmonicsFiveDen[5] = {0, 0, 0, 0, 0};

      if (!fCorrelDenoms[4]) {
        fCorrelDenoms[4] = recursion(5, harmonicsFiveDen).Re();
      }

      TComplex fiveRecursion = recursion(5, harmonicsFiveNum) / fCorrelDenoms[4];

      correlData[0] = fiveRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>
      correlData[1] = fCorrelDenoms[4];   // weight
      correlData[2] = fiveRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>
    } break;
    case 6: {
      int harmonicsSixNum[6] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5]};
      int harmonicsSixDen[6] = {0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[5]) {
        fCorrelDenoms[5] = recursion(6, harmonicsSixDen).Re();
      }

      TComplex sixRecursion = recursion(6, harmonicsSixNum) / fCorrelDenoms[5];

      correlData[0] = sixRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>
      correlData[1] = fCorrelDenoms[5];  // weight
      correlData[2] = sixRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>
    } break;
    case 7: {
      int harmonicsSevenNum[7] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5], harmo[6]};
      int harmonicsSevenDen[7] = {0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[6]) {
        fCorrelDenoms[6] = recursion(7, harmonicsSevenDen).Re();
      }

      TComplex sevenRecursion = recursion(7, harmonicsSevenNum) / fCorrelDenoms[6];

      correlData[0] = sevenRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>
      correlData[1] = fCorrelDenoms[6];    // weight
      correlData[2] = sevenRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>
    } break;
    case 8: {
      int harmonicsEightNum[8] = {harmo[0], harmo[1], harmo[2], harmo[3],
                                  harmo[4], harmo[5], harmo[6], harmo[7]};
      int harmonicsEightDen[8] = {0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[7]) {
        fCorrelDenoms[7] = recursion(8, harmonicsEightDen).Re();
      }

      TComplex eightRecursion = recursion(8, harmonicsEightNum) / fCorrelDenoms[7];

      correlData[0] = eightRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
      correlData[1] = fCorrelDenoms[7];    // weight
      correlData[2] = eightRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
    } break;
    case 9: {
      int harmonicsNineNum[9] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4],
                                 harmo[5], harmo[6], harmo[7], harmo[8]};
      int harmonicsNineDen[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[8]) {
        fCorrelDenoms[8] = recursion(9, harmonicsNineDen).Re();
      }

      TComplex nineRecursion = recursion(9, harmonicsNineNum) / fCorrelDenoms[8];

      correlData[0] = nineRecursion.Re();
      correlData[1] = fCorrelDenoms[8];
      correlData[2] = nineRecursion.Im();
    } break;
    case 10: {
      int harmonicsTenNum[10] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4],
                                 harmo[5], harmo[6], harmo[7], harmo[8], harmo[9]};
      int harmonicsTenDen[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[9]) {
        fCorrelDenoms[9] = recursion(10, harmonicsTenDen).Re();
      }

      TComplex tenRecursion = recursion(10, harmonicsTenNum) / fCorrelDenoms[9];

      correlData[0] = tenRecursion.Re();
      correlData[1] = fCorrelDenoms[9];
      correlData[2] = tenRecursion.Im();
    } break;
    case 12: {
      int harmonicsTwelveNum[12] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5],
                                    harmo[6], harmo[7], harmo[8], harmo[9], harmo[10], harmo[11]};
      int harmonicsTwelveDen[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[11]) {
        fCorrelDenoms[11] = recursion(12, harmonicsTwelveDen).Re();
      }

      TComplex twelveRecursion = recursion(12, harmonicsTwelveNum) / fCorrelDenoms[11];

      correlData[0] = twelveRecursion.Re();
      correlData[1] = fCorrelDenoms[11];
      correlData[2] = twelveRecursion.Im();
    } break;
    case 14: {
      int harmonicsFourteenNum[14] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5], harmo[6],
                                      harmo[7], harmo[8], harmo[9], harmo[10], harmo[11], harmo[12], harmo[13]};
      int harmonicsFourteenDen[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[13]) {
        fCorrelDenoms[13] = recursion(14, harmonicsFourteenDen).Re();
      }

      TComplex fourteenRecursion = recursion(14, harmonicsFourteenNum) / fCorrelDenoms[13];

      correlData[0] = fourteenRecursion.Re();
      correlData[1] = fCorrelDenoms[13];
      correlData[2] = fourteenRecursion.Im();
    } break;
  }
}
int FlowJSPCAnalysis::getCentBin(float cValue)
{
  const float centClasses[] = {0., 1., 2., 5., 10., 20., 30., 40., 50., 60., 70.};
  for (int i = 0; i < 8; i++) {
    if (cValue >= centClasses[i]) {
      continue;
    } else {
      return i - 1;
    }
  }

  // We went through all centrality edges without returning at all.
  // --> The measured percentile is larger than the final class we consider.
  return -1;
}
