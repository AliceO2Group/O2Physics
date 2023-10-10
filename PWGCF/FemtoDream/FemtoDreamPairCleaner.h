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

/// \file FemtoDreamPairCleaner.h
/// \brief FemtoDreamPairCleaner - Makes sure only proper candidates are paired
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de, Laura Serksnyte <laura.serksnyte@cern.ch>, TU München

#ifndef PWGCF_FEMTODREAM_FEMTODREAMPAIRCLEANER_H_
#define PWGCF_FEMTODREAM_FEMTODREAMPAIRCLEANER_H_

#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamPairCleaner
/// \brief Class taking care that no autocorrelations enter the same event distribution
/// \tparam partOne Type of particle 1 (Track/V0/Cascade/...)
/// \tparam partTwo Type of particle 2 (Track/V0/Cascade/...)
template <o2::aod::femtodreamparticle::ParticleType partOne, o2::aod::femtodreamparticle::ParticleType partTwo>
class FemtoDreamPairCleaner
{
 public:
  /// Destructor
  virtual ~FemtoDreamPairCleaner() = default;

  /// Initialization of the QA histograms
  /// \param registry HistogramRegistry
  void init(HistogramRegistry* registry)
  {
    if (registry) {
      mHistogramRegistry = registry;
      // \todo some QA histograms like in FemtoDream
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
    if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
      /// Track-Track combination
      if (part1.partType() != o2::aod::femtodreamparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtodreamparticle::ParticleType::kTrack) {
        LOG(fatal) << "FemtoDreamPairCleaner: passed arguments don't agree with FemtoDreamPairCleaner instantiation! Please provide kTrack,kTrack candidates.";
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kV0) {
      /// Track-V0 combination
      if (part2.partType() != o2::aod::femtodreamparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoDreamPairCleaner: passed arguments don't agree with FemtoDreamPairCleaner instantiation! Please provide second argument kV0 candidate.";
        return false;
      }
      const auto& posChild = particles.iteratorAt(part2.index() - 2);
      const auto& negChild = particles.iteratorAt(part2.index() - 1);
      if (part1.globalIndex() != posChild.globalIndex() || part2.globalIndex() != negChild.globalIndex()) {
        return true;
      }
      return false;
    } else {
      LOG(fatal) << "FemtoDreamPairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                             ///< For QA output
  static constexpr o2::aod::femtodreamparticle::ParticleType mPartOneType = partOne; ///< Type of particle 1
  static constexpr o2::aod::femtodreamparticle::ParticleType mPartTwoType = partTwo; ///< Type of particle 2
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_FEMTODREAMPAIRCLEANER_H_
