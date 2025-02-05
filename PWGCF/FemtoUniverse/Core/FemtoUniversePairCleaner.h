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

/// \file FemtoUniversePairCleaner.h
/// \brief FemtoUniversePairCleaner - Makes sure only proper candidates are paired
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Laura Serksnyte, TU München,laura.serksnyte@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRCLEANER_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRCLEANER_H_

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

namespace o2::analysis::femto_universe
{

/// \class FemtoUniversePairCleaner
/// \brief Class taking care that no autocorrelations enter the same event distribution
/// \tparam partOne Type of particle 1 (Track/V0/Cascade/...)
/// \tparam partTwo Type of particle 2 (Track/V0/Cascade/...)
template <o2::aod::femtouniverseparticle::ParticleType partOne, o2::aod::femtouniverseparticle::ParticleType partTwo>
class FemtoUniversePairCleaner
{
 public:
  /// Destructor
  virtual ~FemtoUniversePairCleaner() = default;

  /// Initialization of the QA histograms
  /// \param registry HistogramRegistry
  void init(HistogramRegistry* registry)
  {
    if (registry) {
      mHistogramRegistry = registry;
      // \todo some QA histograms like in FemtoUniverse
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
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      /// Track-Track combination
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide kTrack,kTrack candidates.";
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack) {
      /// Track-Track combination
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide kMCTruthTrack,kMCTruthTrack candidates.";
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      /// Track-V0 combination part1 is hadron and part2 is v0
      if (part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide second argument kV0 candidate.";
        return false;
      }
      // Getting v0 (part2) children
      const auto& posChild = particles.iteratorAt(part2.index() - 2);
      const auto& negChild = particles.iteratorAt(part2.index() - 1);
      if (part1.globalIndex() == posChild.globalIndex() || part1.globalIndex() == negChild.globalIndex()) {
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kV0 && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      /// V0-V0 combination both part1 and part2 are v0
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0 || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide first and second arguments kV0 candidate.";
        return false;
      }
      // Getting v0 children for part1
      const auto& posChild1 = particles.iteratorAt(part1.index() - 2);
      const auto& negChild1 = particles.iteratorAt(part1.index() - 1);
      // Getting v0 children for part2
      const auto& posChild2 = particles.iteratorAt(part2.index() - 2);
      const auto& negChild2 = particles.iteratorAt(part2.index() - 1);
      if (posChild1.globalIndex() == posChild2.globalIndex() || negChild1.globalIndex() == negChild2.globalIndex()) {
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// Track-Cascade combination part1 is hadron and part2 is cascade
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kCascade) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide first argument kTrack candidate and second argument kCascade candidate.";
        return false;
      }
      // Getting cascade children for part2
      const auto& posChild = particles.iteratorAt(part2.index() - 3);
      const auto& negChild = particles.iteratorAt(part2.index() - 2);
      const auto& bachelor = particles.iteratorAt(part2.index() - 1);
      if (part1.globalIndex() == posChild.globalIndex() || part1.globalIndex() == negChild.globalIndex() || part1.globalIndex() == bachelor.globalIndex()) {
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kCascade && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// Cascade-Cascade combination both part1 and part2 are cascades
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kCascade || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kCascade) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide first and second arguments kCascade candidate.";
        return false;
      }
      // Getting cascade children for part1
      const auto& posChild1 = particles.iteratorAt(part1.index() - 3);
      const auto& negChild1 = particles.iteratorAt(part1.index() - 2);
      const auto& bachelor1 = particles.iteratorAt(part1.index() - 1);
      // Getting cascade children for part2
      const auto& posChild2 = particles.iteratorAt(part2.index() - 3);
      const auto& negChild2 = particles.iteratorAt(part2.index() - 2);
      const auto& bachelor2 = particles.iteratorAt(part2.index() - 1);
      if (posChild1.globalIndex() == posChild2.globalIndex() || negChild1.globalIndex() == negChild2.globalIndex() || bachelor1.globalIndex() == bachelor2.globalIndex()) {
        return false;
      }
      return part1.globalIndex() != part2.globalIndex();
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kD0) {
      /// Track-D0 combination part1 is hadron and part2 is D0
      if (part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kD0) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide second argument kD0 candidate.";
        return false;
      }
      // Getting D0 (part2) children
      const auto& posChild = particles.iteratorAt(part2.index() - 2);
      const auto& negChild = particles.iteratorAt(part2.index() - 1);

      if (part1.globalIndex() != posChild.globalIndex() && part1.globalIndex() != negChild.globalIndex()) {
        return true;
      }
      return false;
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
      /// Track-Phi combination part1 is Phi and part 2 is hadron
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kPhi) {
        LOG(fatal) << "FemtoUniversePairCleaner: passed arguments don't agree with FemtoUniversePairCleaner instantiation! Please provide second argument kPhi candidate.";
        return false;
      }

      // getting Phi (part1) children
      const auto& posChild = particles.iteratorAt(part2.index() - 2);
      const auto& negChild = particles.iteratorAt(part2.index() - 1);

      if (part1.globalIndex() != posChild.globalIndex() && part1.globalIndex() != negChild.globalIndex()) {
        return true;
      }
      return false;
    } else {
      LOG(fatal) << "FemtoUniversePairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                ///< For QA output
  static constexpr o2::aod::femtouniverseparticle::ParticleType kPartOneType = partOne; ///< Type of particle 1
  static constexpr o2::aod::femtouniverseparticle::ParticleType kPartTwoType = partTwo; ///< Type of particle 2
};
} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPAIRCLEANER_H_
