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

/// \file ParticleOrigin.h
/// \brief commonly used enums for particle origins.
/// \author marvin.hemmer@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_PARTICLEORIGIN_H_
#define PWGEM_PHOTONMESON_UTILS_PARTICLEORIGIN_H_

#include <cstdint>

namespace o2::analysis::em
{

// Classifies the production history of a lepton
enum class LeptonOrigin : uint8_t {
  Conversion = 0,
  Compton = 1,
  PhotoElectric = 2,
  DeltaRay = 3,
  DirectMesonDecay = 4,
  Other = 5
};

// Classifies the production history of a photon
enum class PhotonOrigin : uint8_t {
  Direct = 0,         // direct photons from quarks or gluons from inital scattering
  Decay = 1,          // decay photons from
  Bremsstrahlung = 2, // photon from bremsstrahlung
  Annihilation = 3,   // e+e- -> gammagamma
  Hadronic = 4,       // hadronic interactions like charged pions with material
  Other = 5           // anything else
};

} // namespace o2::analysis::em

#endif // PWGEM_PHOTONMESON_UTILS_PARTICLEORIGIN_H_
