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

/// \file FemtoWorldPairCleaner.h
/// \brief FemtoWorldPairCleaner - Makes sure only proper candidates are paired
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de, Laura Serksnyte <laura.serksnyte@cern.ch>, TU München
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch

#ifndef PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPAIRCLEANER_H_
#define PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPAIRCLEANER_H_

#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoWorld
{

/// \class FemtoWorldPairCleaner
/// \brief Class taking care that no autocorrelations enter the same event distribution
/// \tparam partOne Type of particle 1 (Track/V0/Cascade/...)
/// \tparam partTwo Type of particle 2 (Track/V0/Cascade/...)
template <o2::aod::femtoworldparticle::ParticleType partOne, o2::aod::femtoworldparticle::ParticleType partTwo>
class FemtoWorldPairCleaner
{
 public:
  /// Destructor
  virtual ~FemtoWorldPairCleaner() = default;

  /// Initialization of the QA histograms
  /// \param registry HistogramRegistry
  void init(HistogramRegistry* registry)
  {
    if (registry) {
      mHistogramRegistry = registry;
      // \todo some QA histograms like in FemtoWorld
    }
  }

  /// Check whether a given pair has shared tracks
  /// \tparam Part Data type of the particle
  /// \tparam Parts Data type of the collection of all particles
  /// \param part1 Particle 1
  /// \param part2 Particle 2
  /// \param particles Collection of all particles passed to the task
  /// \return Whether the pair has shared tracks
  template <typename Part, typename Parts>
  bool isCleanPair(Part const& part1, Part const& part2, Parts const& particles)
  {
    if constexpr (mPartOneType == o2::aod::femtoworldparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtoworldparticle::ParticleType::kTrack) {
      /// Track-Track combination
      if (part1.partType() != o2::aod::femtoworldparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtoworldparticle::ParticleType::kTrack) {
        LOG(fatal) << "FemtoWorldPairCleaner: passed arguments don't agree with FemtoWorldPairCleaner instantiation! Please provide kTrack,kTrack candidates.";
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (mPartOneType == o2::aod::femtoworldparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtoworldparticle::ParticleType::kV0) {
      /// Track-V0 combination
      if (part2.partType() != o2::aod::femtoworldparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoWorldPairCleaner: passed arguments don't agree with FemtoWorldPairCleaner instantiation! Please provide second argument kV0 candidate.";
        return false;
      }
      uint64_t id1 = part2.index() - 2;
      uint64_t id2 = part2.index() - 1;
      auto daughter1 = particles.begin() + id1;
      auto daughter2 = particles.begin() + id2;
      if ((*daughter1).indices()[0] <= 0 && (*daughter1).indices()[1] <= 0 && (*daughter2).indices()[0] <= 0 && (*daughter2).indices()[1] <= 0) {
        return true;
      }
      return false;
    } else if constexpr (mPartOneType == o2::aod::femtoworldparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtoworldparticle::ParticleType::kPhi) {
      return true;
    } else if constexpr (mPartOneType == o2::aod::femtoworldparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtoworldparticle::ParticleType::kD0D0bar) {
      return true;
    } else {
      LOG(fatal) << "FemtoWorldPairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                             ///< For QA output
  static constexpr o2::aod::femtoworldparticle::ParticleType mPartOneType = partOne; ///< Type of particle 1
  static constexpr o2::aod::femtoworldparticle::ParticleType mPartTwoType = partTwo; ///< Type of particle 2
};
} // namespace o2::analysis::femtoWorld

#endif // PWGCF_FEMTOWORLD_CORE_FEMTOWORLDPAIRCLEANER_H_
