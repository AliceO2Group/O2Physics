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

/// \file FemtoUtils.h
/// \brief Utilities for the FemtoDream framework
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#ifndef PWGCF_FEMTODREAM_FEMTOUTILS_H_
#define PWGCF_FEMTODREAM_FEMTOUTILS_H_

#include <vector>
#include <functional>
#include <algorithm>
#include "Framework/ASoAHelpers.h"
#include "PWGCF/DataModel/FemtoDerived.h"

namespace o2::analysis::femtoDream
{

enum kDetector { kTPC,
                 kTPCTOF,
                 kNdetectors };

/// internal function that returns the kPIDselection element corresponding to a
/// specifica n-sigma value \param nSigma number of sigmas for PID
/// \param vNsigma vector with the number of sigmas of interest
/// \return kPIDselection corresponding to n-sigma
int getPIDselection(float nSigma, std::vector<float> vNsigma)
{
  std::sort(vNsigma.begin(), vNsigma.end(), std::greater<>());
  auto it = std::find(vNsigma.begin(), vNsigma.end(), nSigma);
  if (it == vNsigma.end()) {
    it = vNsigma.begin() + 1;
    LOG(warn) << "Invalid value of nSigma: " << nSigma << ". Return the first value of the vector: " << *(it);
  }
  return std::distance(vNsigma.begin(), it);
}

/// function that checks whether the PID selection specified in the vectors is
/// fulfilled
/// \param pidcut Bit-wise container for the PID
/// \param vSpecies vector with ID corresponding to the selected species (output from cutculator)
/// \param nSpecies number of available selected species (output from cutculator)
/// \param nSigma number of sigma selection fo PID
/// \param vNsigma vector with available n-sigma selections for PID
/// \param kDetector enum corresponding to the PID technique
/// \return Whether the PID selection specified in the vectors is fulfilled
bool isPIDSelected(aod::femtodreamparticle::cutContainerType pidcut,
                   int vSpecies,
                   int nSpecies,
                   float nSigma,
                   std::vector<float> vNsigma,
                   kDetector iDet)
{
  int iNsigma = getPIDselection(nSigma, vNsigma);
  int nDet = static_cast<int>(kDetector::kNdetectors);
  int bit_to_check = 1 + (vNsigma.size() - (iNsigma + 1)) * nDet * nSpecies + (nSpecies - (vSpecies + 1)) * nSpecies + (nDet - 1 - iDet);
  return ((pidcut >> (bit_to_check)) & 1) == 1;
};

/// function that checks whether the PID selection specified in the vectors is fulfilled, depending on the momentum TPC or TPC+TOF PID is conducted
/// \param pidcut Bit-wise container for the PID
/// \param momentum Momentum of the track
/// \param pidThresh Momentum threshold that separates between TPC and TPC+TOF PID
/// \param vSpecies Vector with the species of interest (number returned by the CutCulator)
/// \param nSpecies number of available selected species (output from cutculator)
/// \param nSigmaTPC Number of TPC sigmas for selection
/// \param nSigmaTPCTOF Number of TPC+TOF sigmas for selection (circular selection)
/// \return Whether the PID selection is fulfilled
bool isFullPIDSelected(aod::femtodreamparticle::cutContainerType const& pidCut,
                       float momentum,
                       float pidThresh,
                       int vSpecies,
                       int nSpecies,
                       std::vector<float> vNsigma,
                       float nSigmaTPC,
                       float nSigmaTPCTOF)
{
  bool pidSelection = true;
  if (momentum < pidThresh) {
    /// TPC PID only
    pidSelection = isPIDSelected(pidCut, vSpecies, nSpecies, nSigmaTPC, vNsigma, kDetector::kTPC);
  } else {
    /// TPC + TOF PID
    pidSelection = isPIDSelected(pidCut, vSpecies, nSpecies, nSigmaTPCTOF, vNsigma, kDetector::kTPCTOF);
  }
  return pidSelection;
};

int checkDaughterType(o2::aod::femtodreamparticle::ParticleType partType, int motherPDG)
{
  int partOrigin = 0;
  if (partType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
    switch (abs(motherPDG)) {
      case 3122:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterLambda;
        break;
      case 3222:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterSigmaplus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    } // switch

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kV0) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kV0Child) {
    switch (abs(motherPDG)) {
      case 3122:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterLambda;
        break;
      case 3222:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterSigmaplus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    } // switch

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
  }
  return partOrigin;
};

} // namespace o2::analysis::femtoDream
#endif // PWGCF_FEMTODREAM_FEMTOUTILS_H_
