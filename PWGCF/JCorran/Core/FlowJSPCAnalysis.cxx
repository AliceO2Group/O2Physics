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

TComplex FlowJSPCAnalysis::Q(const Int_t harmN, const Int_t p)
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
TComplex FlowJSPCAnalysis::Recursion(int n, int* harmonic, int mult = 1, int skip = 0)
{
  Int_t nm1 = n - 1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0)
    return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip)
    return c;

  Int_t multp1 = mult + 1;
  Int_t nm2 = n - 2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n - 3;

  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1)
    return c - c2;
  return c - Double_t(mult) * c2;
} // End of recursion

void FlowJSPCAnalysis::CalculateCorrelators(const Int_t fCentBin)
{
  // Loop over the combinations of harmonics and calculate the corresponding SPC num and den.

  // Declare the arrays to later fill all the needed bins for the correlators
  // and the error terms.
  Double_t* dataCorrelation = new Double_t[3]; // cosine, weight, sine.
  Double_t correlationNum;
  Double_t weightCorrelationNum;
  Double_t correlationDenom;
  Double_t weightCorrelationDenom;

  for (Int_t j = 0; j < 12; j++) {
    if (fHarmosArray[j][0] == 0) {
      continue;
    } // Skip null correlator list.

    // Calculate the numerator.
    Int_t hArrayNum[7] = {0};
    for (int iH = 0; iH < 7; iH++) {
      hArrayNum[iH] = fHarmosArray[j][iH + 1];
    }
    Correlation(fHarmosArray[j][0], 7, hArrayNum, dataCorrelation);
    correlationNum = dataCorrelation[0];
    weightCorrelationNum = dataCorrelation[1];

    // Calculate the denominator.
    Int_t nPartDen = 2 * fHarmosArray[j][0];
    Int_t hArrayDen[14] = {0};
    for (int iH = 0; iH < 7; iH++) {
      hArrayDen[2 * iH] = hArrayNum[iH];
      hArrayDen[2 * iH + 1] = -1 * hArrayNum[iH];
    }

    Correlation(nPartDen, 14, hArrayDen, dataCorrelation);
    correlationDenom = dataCorrelation[0];
    weightCorrelationDenom = dataCorrelation[1];

    // Histogram filling

    FillHistograms(fCentBin, j, correlationNum, correlationDenom, weightCorrelationNum, weightCorrelationDenom);

    correlationNum = 0.;
    weightCorrelationNum = 0.;
    correlationDenom = 0.;
    weightCorrelationDenom = 0.;
  }
}

void FlowJSPCAnalysis::FillHistograms(const Int_t fCentBin, Int_t ind, Double_t cNum, Double_t cDenom, Double_t wNum, Double_t wDenom)
{
  switch (fCentBin) {
    case 0: {
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    case 1: {
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    case 2: {
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    case 3: {
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    case 4: {
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    case 5: {
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    case 6: {
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    case 7: {
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fResults"), 2. * (Float_t)(ind) + 0.5, cNum, wNum);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fResults"), 2. * (Float_t)(ind) + 1.5, cDenom, wDenom);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 0.5, cNum * cDenom, wNum * wDenom);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 1.5, wNum * wDenom, 1.);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 2.5, wNum, 1.);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fCovResults"), 2. * (Float_t)(ind) + 3.5, wDenom, 1.);
    } break;
    default:
      return;
  }
}

void FlowJSPCAnalysis::FillQAHistograms(const Int_t fCentBin, Double_t phi, Double_t phiWeight)
{
  switch (fCentBin) {
    case 0: {
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 1: {
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 2: {
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 3: {
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 4: {
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 5: {
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 6: {
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    case 7: {
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("phiBefore"), phi);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("phiAfter"), phi, phiWeight);
    } break;
    default:
      return;
  }
}

void FlowJSPCAnalysis::Correlation(Int_t c_nPart, Int_t c_nHarmo, Int_t* harmo, Double_t* correlData)
{
  // Calculate the correlators for the provided set of harmonics using Q-vectors.
  // Protection against anisotropic correlators.
  Int_t sumHarmo = 0;
  for (Int_t i = 0; i < c_nHarmo; i++) {
    sumHarmo += harmo[i];
  }
  if (sumHarmo != 0) {
    printf("\nOups, this correlator is not isotropic(sum = %d). Bye\n", sumHarmo);
    return;
  }

  switch (c_nPart) {
    case 2: {
      Int_t harmonicsTwoNum[2] = {harmo[0], harmo[1]};
      Int_t harmonicsTwoDen[2] = {0, 0};

      if (!fCorrelDenoms[1]) {
        fCorrelDenoms[1] = Recursion(2, harmonicsTwoDen).Re();
      }

      TComplex twoRecursion = Recursion(2, harmonicsTwoNum) / fCorrelDenoms[1];

      correlData[0] = twoRecursion.Re(); // <cos(h1*phi1+h2*phi2)>
      correlData[1] = fCorrelDenoms[1];  // weight
      correlData[2] = twoRecursion.Im(); // <sin(h1*phi1+h2*phi2)>
    } break;
    case 3: {
      Int_t harmonicsThreeNum[3] = {harmo[0], harmo[1], harmo[2]};
      Int_t harmonicsThreeDen[3] = {0, 0, 0};

      if (!fCorrelDenoms[2]) {
        fCorrelDenoms[2] = Recursion(3, harmonicsThreeDen).Re();
      }

      TComplex threeRecursion = Recursion(3, harmonicsThreeNum) / fCorrelDenoms[2];

      correlData[0] = threeRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3)>
      correlData[1] = fCorrelDenoms[2];    // weight
      correlData[2] = threeRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3)>
    } break;
    case 4: {
      Int_t harmonicsFourNum[4] = {harmo[0], harmo[1], harmo[2], harmo[3]};
      Int_t harmonicsFourDen[4] = {0, 0, 0, 0};

      if (!fCorrelDenoms[3]) {
        fCorrelDenoms[3] = Recursion(4, harmonicsFourDen).Re();
      }

      TComplex fourRecursion = Recursion(4, harmonicsFourNum) / fCorrelDenoms[3];

      correlData[0] = fourRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>
      correlData[1] = fCorrelDenoms[3];   // weight
      correlData[2] = fourRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>
    } break;
    case 5: {
      Int_t harmonicsFiveNum[5] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4]};
      Int_t harmonicsFiveDen[5] = {0, 0, 0, 0, 0};

      if (!fCorrelDenoms[4]) {
        fCorrelDenoms[4] = Recursion(5, harmonicsFiveDen).Re();
      }

      TComplex fiveRecursion = Recursion(5, harmonicsFiveNum) / fCorrelDenoms[4];

      correlData[0] = fiveRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>
      correlData[1] = fCorrelDenoms[4];   // weight
      correlData[2] = fiveRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>
    } break;
    case 6: {
      Int_t harmonicsSixNum[6] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5]};
      Int_t harmonicsSixDen[6] = {0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[5]) {
        fCorrelDenoms[5] = Recursion(6, harmonicsSixDen).Re();
      }

      TComplex sixRecursion = Recursion(6, harmonicsSixNum) / fCorrelDenoms[5];

      correlData[0] = sixRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>
      correlData[1] = fCorrelDenoms[5];  // weight
      correlData[2] = sixRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>
    } break;
    case 7: {
      Int_t harmonicsSevenNum[7] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5], harmo[6]};
      Int_t harmonicsSevenDen[7] = {0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[6]) {
        fCorrelDenoms[6] = Recursion(7, harmonicsSevenDen).Re();
      }

      TComplex sevenRecursion = Recursion(7, harmonicsSevenNum) / fCorrelDenoms[6];

      correlData[0] = sevenRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>
      correlData[1] = fCorrelDenoms[6];    // weight
      correlData[2] = sevenRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>
    } break;
    case 8: {
      Int_t harmonicsEightNum[8] = {harmo[0], harmo[1], harmo[2], harmo[3],
                                    harmo[4], harmo[5], harmo[6], harmo[7]};
      Int_t harmonicsEightDen[8] = {0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[7]) {
        fCorrelDenoms[7] = Recursion(8, harmonicsEightDen).Re();
      }

      TComplex eightRecursion = Recursion(8, harmonicsEightNum) / fCorrelDenoms[7];

      correlData[0] = eightRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
      correlData[1] = fCorrelDenoms[7];    // weight
      correlData[2] = eightRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
    } break;
    case 9: {
      Int_t harmonicsNineNum[9] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4],
                                   harmo[5], harmo[6], harmo[7], harmo[8]};
      Int_t harmonicsNineDen[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[8]) {
        fCorrelDenoms[8] = Recursion(9, harmonicsNineDen).Re();
      }

      TComplex nineRecursion = Recursion(9, harmonicsNineNum) / fCorrelDenoms[8];

      correlData[0] = nineRecursion.Re();
      correlData[1] = fCorrelDenoms[8];
      correlData[2] = nineRecursion.Im();
    } break;
    case 10: {
      Int_t harmonicsTenNum[10] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4],
                                   harmo[5], harmo[6], harmo[7], harmo[8], harmo[9]};
      Int_t harmonicsTenDen[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[9]) {
        fCorrelDenoms[9] = Recursion(10, harmonicsTenDen).Re();
      }

      TComplex tenRecursion = Recursion(10, harmonicsTenNum) / fCorrelDenoms[9];

      correlData[0] = tenRecursion.Re();
      correlData[1] = fCorrelDenoms[9];
      correlData[2] = tenRecursion.Im();
    } break;
    case 12: {
      Int_t harmonicsTwelveNum[12] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5],
                                      harmo[6], harmo[7], harmo[8], harmo[9], harmo[10], harmo[11]};
      Int_t harmonicsTwelveDen[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[11]) {
        fCorrelDenoms[11] = Recursion(12, harmonicsTwelveDen).Re();
      }

      TComplex twelveRecursion = Recursion(12, harmonicsTwelveNum) / fCorrelDenoms[11];

      correlData[0] = twelveRecursion.Re();
      correlData[1] = fCorrelDenoms[11];
      correlData[2] = twelveRecursion.Im();
    } break;
    case 14: {
      Int_t harmonicsFourteenNum[14] = {harmo[0], harmo[1], harmo[2], harmo[3], harmo[4], harmo[5], harmo[6],
                                        harmo[7], harmo[8], harmo[9], harmo[10], harmo[11], harmo[12], harmo[13]};
      Int_t harmonicsFourteenDen[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      if (!fCorrelDenoms[13]) {
        fCorrelDenoms[13] = Recursion(14, harmonicsFourteenDen).Re();
      }

      TComplex fourteenRecursion = Recursion(14, harmonicsFourteenNum) / fCorrelDenoms[13];

      correlData[0] = fourteenRecursion.Re();
      correlData[1] = fCorrelDenoms[13];
      correlData[2] = fourteenRecursion.Im();
    } break;
  }
}
int FlowJSPCAnalysis::GetCentBin(float cValue)
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
