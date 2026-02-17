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

/// \file EMNonLin.h
/// \brief Header file of NonLin class for photons.
/// \author M. Hemmer, marvin.hemmer@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_EMNONLIN_H_
#define PWGEM_PHOTONMESON_CORE_EMNONLIN_H_

#include <Framework/Configurable.h>

#include <cstdint> // uint8_t

namespace o2::pwgem::nonlin
{

/// \class EMNonLin
/// \brief Dynamically-sized bit container with bit-level storage.
///
/// Bits can be set beyond the current size, in which case the container
/// grows automatically. Access via test() and reset() requires the index
/// to be within the current size. Bits are all on by default and will be
/// switched off if particle fails a cut
class EMNonLin
{
 public:
  enum class PhotonType : uint8_t {
    kEMC = 0,
    kPCM = 1,
    kPHOS = 2 // just in case
  };

  /// \brief gets the correction value for energy or pT for a specicif
  /// \param inputCalibValue pT or energy of the photon that needs calibration
  /// \param photonType type of the photon (e.g. 0 for EMC)
  /// \param cent centrailty of the current collision in case the correction is centrality dependent
  float getCorrectionFactor(float inputCalibValue, PhotonType photonType = PhotonType::kEMC, float cent = 0);

 private:
};
} // namespace o2::pwgem::nonlin

#endif // PWGEM_PHOTONMESON_CORE_EMNONLIN_H_
