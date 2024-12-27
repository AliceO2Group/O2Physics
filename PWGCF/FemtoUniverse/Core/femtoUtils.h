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

/// \file femtoUtils.h
/// \brief Utilities for the FemtoUniverse framework
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUTILS_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUTILS_H_

#include <vector>
#include <functional>
#include <algorithm>
#include "Framework/ASoAHelpers.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe
{

enum KDetector { kTPC,
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
    LOG(warn) << "Invalid value of nSigma: " << nSigma << ". Return the largest value of the vector: " << *it;
  }
  return std::distance(vNsigma.begin(), it);
}

/// function that checks whether the PID selection specified in the vectors is
/// fulfilled
/// \param pidcut Bit-wise container for the PID
/// \param vSpecies number (ID) of the species of interest (as returned by the CutCulator), e.g. 0 / 1 / 2, usually 0 if only one particle species in the skimmed data
/// \param nSpecies number of available selected species (output from cutculator), i.e. how many particle types were saved in the skimmed data
/// \param nSigma Nsigma selection for PID (e.g. 3, for NsigmaTPC < 3 or NsigmaTPCTOF < 3)
/// \param vNsigma vector with available n-sigma selections for PID (to check if chosen nSigma value is avialable + size to get the bit number)
/// \param KDetector enum corresponding to the PID technique
/// \return Whether the PID selection specified in the vectors is fulfilled
bool isPIDSelected(aod::femtouniverseparticle::CutContainerType pidcut,
                   int vSpecies,
                   int nSpecies,
                   float nSigma,
                   std::vector<float> vNsigma,
                   KDetector iDet)
{
  int iNsigma = getPIDselection(nSigma, vNsigma);
  int nDet = static_cast<int>(KDetector::kNdetectors);
  int bitToCheck = 1 + (vNsigma.size() - (iNsigma + 1)) * nDet * nSpecies + (nSpecies - (vSpecies + 1)) * nSpecies + (nDet - 1 - iDet);
  return ((pidcut >> (bitToCheck)) & 1) == 1;
};

/// function that checks whether the PID selection specified in the vectors is fulfilled, depending on the momentum TPC or TPC+TOF PID is conducted
/// \param pidcut Bit-wise container for the PID
/// \param momentum Momentum of the track
/// \param pidThresh Momentum threshold that separates between TPC and TPC+TOF PID
/// \param vSpecies number (ID) of the species of interest (as returned by the CutCulator), e.g. 0 / 1 / 2, usually 0 if only one particle specie in the skimmed data
/// \param nSpecies number of available selected species (output from cutculator), i.e. how many particle types were saved in the skimmed data
/// \param vNsigma vector with available n-sigma selections for PID
/// \param nSigmaTPC Number of TPC sigmas for selection
/// \param nSigmaTPCTOF Number of TPC+TOF sigmas for selection (circular selection)
/// \return Whether the PID selection is fulfilled
bool isFullPIDSelected(aod::femtouniverseparticle::CutContainerType const& pidCut,
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
    pidSelection = isPIDSelected(pidCut, vSpecies, nSpecies, nSigmaTPC, vNsigma, KDetector::kTPC);
  } else {
    /// TPC + TOF PID
    pidSelection = isPIDSelected(pidCut, vSpecies, nSpecies, nSigmaTPCTOF, vNsigma, KDetector::kTPCTOF);
  }
  return pidSelection;
};

int checkDaughterType(o2::aod::femtouniverseparticle::ParticleType partType, int motherPDG)
{
  int partOrigin = 0;
  if (partType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {

    switch (std::abs(motherPDG)) {
      case 3122:
        partOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kDaughterLambda;
        break;
      case 3222:
        partOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kDaughterSigmaplus;
        break;
      default:
        partOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kDaughter;
    } // switch

  } else if (partType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
    partOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kDaughter;

  } else if (partType == o2::aod::femtouniverseparticle::ParticleType::kV0Child) {
    partOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kDaughter;

  } else if (partType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
    partOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kDaughter;

  } else if (partType == o2::aod::femtouniverseparticle::ParticleType::kCascadeBachelor) {
    partOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kDaughter;
  }
  return partOrigin;
};

} // namespace o2::analysis::femto_universe
#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUTILS_H_
