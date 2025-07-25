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

/// \file CFFilter.cxx
/// \brief Utilities for the FemtoWorld framework
///
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch

#ifndef FEMTOWORLD_UTILS_H_
#define FEMTOWORLD_UTILS_H_

#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"

#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>

#include <vector>

namespace o2::analysis::femtoWorld
{

enum kPIDselection {
  k3d5sigma = 0,
  k3sigma = 1,
  k2d5sigma = 2
};

enum kDetector {
  kTPC = 0,
  kTPCTOF = 1,
  kNdetectors = 2
};

/// internal function that returns the kPIDselection element corresponding to a specifica n-sigma value
/// \param nSigma number of sigmas for PID
/// \param vNsigma vector with the number of sigmas of interest
/// \return kPIDselection corresponding to n-sigma
kPIDselection getPIDselection(const float nSigma, const std::vector<float>& vNsigma)
{
  for (int i = 0; i < (int)vNsigma.size(); i++) {
    if (abs(nSigma - vNsigma[i]) < 1e-3) {
      return static_cast<kPIDselection>(i);
    }
  }
  LOG(info) << "Invalid value of nSigma: " << nSigma << ". Standard 3 sigma returned." << std::endl;
  return kPIDselection::k3sigma;
}

/// function that checks whether the PID selection specified in the vectors is fulfilled
/// \param pidcut Bit-wise container for the PID
/// \param vSpecies vector with ID corresponding to the selected species (output from cutculator)
/// \param nSpecies number of available selected species (output from cutculator)
/// \param nSigma number of sigma selection fo PID
/// \param vNsigma vector with available n-sigma selections for PID
/// \param kDetector enum corresponding to the PID technique
/// \return Whether the PID selection specified in the vectors is fulfilled
bool isPIDSelected(aod::femtoworldparticle::cutContainerType const& pidcut, std::vector<int> const& vSpecies, int nSpecies, float nSigma, const std::vector<float>& vNsigma, const kDetector iDet = kDetector::kTPC)
{
  bool pidSelection = true;
  kPIDselection iNsigma = getPIDselection(nSigma, vNsigma);
  for (auto iSpecies : vSpecies) {
    //\todo we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    // if (!((pidcut >> it.first) & it.second)) {
    int bit_to_check = nSpecies * kDetector::kNdetectors * iNsigma + iSpecies * kDetector::kNdetectors + iDet;
    if (!(pidcut & (1UL << bit_to_check))) {
      pidSelection = false;
    }
  }
  return pidSelection;
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
bool isFullPIDSelected(aod::femtoworldparticle::cutContainerType const& pidCut, float const momentum, float const pidThresh, std::vector<int> const& vSpecies, int nSpecies, const std::vector<float>& vNsigma = {3.5, 3., 2.5}, const float nSigmaTPC = 3.5, const float nSigmaTPCTOF = 3.5)
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

} // namespace o2::analysis::femtoWorld
#endif
