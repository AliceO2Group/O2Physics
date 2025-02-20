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

/// \file FemtoUniverseSoftPionRemoval.h
/// \brief FemtoUniverseSoftPionRemoval - Checking the soft pions from D* decay and removing them
/// \author Katarzyna Gwiździel, WUT, katarzyna.gwizdziel@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESOFTPIONREMOVAL_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESOFTPIONREMOVAL_H_

#include <vector>

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

namespace o2::analysis::femto_universe
{

/// \class FemtoUniverseSoftPionRemoval
/// \brief Class taking care of removing soft pions from D* decays
/// \tparam partOne Type of particle 1 (Track/D0/...)
/// \tparam partTwo Type of particle 2 (Track/D0/...)
template <o2::aod::femtouniverseparticle::ParticleType partOne, o2::aod::femtouniverseparticle::ParticleType partTwo>
class FemtoUniverseSoftPionRemoval
{
 public:
  /// Destructor
  virtual ~FemtoUniverseSoftPionRemoval() = default;

  /// Initialization of the QA histograms
  /// \param registry HistogramRegistry
  void init(HistogramRegistry* registry)
  {
    if (registry) {
      mHistogramRegistry = registry;
      mHistogramRegistry->add("SoftPion/softPionMassVsPt", "; M(K#pi#pi-K#pi); p_{T}", kTH2F, {{200, 0.0, 0.2}, {36, 0., 36.}});
    }
  }

  /// Check whether a track is a soft pion from D* decay
  /// \tparam Part Data type of the particle
  /// \tparam Parts Data type of the collection of all particles
  /// \param part1 Particle 1
  /// \param part2 Particle 2
  /// \param particles Collection of all particles passed to the task
  /// \return Whether the track is a soft pion
  template <typename Part, typename Parts>
  bool isSoftPion(Part const& part1, Part const& part2, Parts const& particles, bool isD0Cand, bool isD0barCand, double sigma)
  {
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kD0) {
      /// Track-D0 combination part1 is hadron and part2 is D0
      if (part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kD0) {
        LOG(fatal) << "FemtoUniverseSoftPionRemoval: passed arguments don't agree with FemtoUniverseSoftPionRemoval instantiation! Please provide second argument kD0 candidate.";
        return false;
      }
      // Getting D0 (part2) children
      const auto& posChild = particles.iteratorAt(part2.index() - 2);
      const auto& negChild = particles.iteratorAt(part2.index() - 1);
      // Pion and kaon mass
      double massPion = o2::constants::physics::MassPiPlus;
      double massKaon = o2::constants::physics::MassKPlus;
      // D* reconstruction
      double pSum2 = std::pow(posChild.px() + negChild.px() + part1.px(), 2.0) + std::pow(posChild.py() + negChild.py() + part1.py(), 2.0) + std::pow(posChild.pz() + negChild.pz() + part1.pz(), 2.0);
      // Energies of the daughters -> D0->K-pi+
      double e1Pi = std::sqrt(std::pow(massPion, 2.0) + std::pow(posChild.px(), 2.0) + std::pow(posChild.py(), 2.0) + std::pow(posChild.pz(), 2.0));
      double e1K = std::sqrt(std::pow(massKaon, 2.0) + std::pow(negChild.px(), 2.0) + std::pow(negChild.py(), 2.0) + std::pow(negChild.pz(), 2.0));
      // Energies of the daughters -> D0bar->K+pi-
      double e2Pi = std::sqrt(std::pow(massPion, 2.0) + std::pow(negChild.px(), 2.0) + std::pow(negChild.py(), 2.0) + std::pow(negChild.pz(), 2.0));
      double e2K = std::sqrt(std::pow(massKaon, 2.0) + std::pow(posChild.px(), 2.0) + std::pow(posChild.py(), 2.0) + std::pow(posChild.pz(), 2.0));
      // Soft pion energy
      auto ePion = RecoDecay::e(massPion, part1.p());
      // D* masses
      double mDstar1 = std::sqrt(std::pow(e1Pi + e1K + ePion, 2.0) - pSum2);
      double mDstar2 = std::sqrt(std::pow(e2Pi + e2K + ePion, 2.0) - pSum2);

      bool isSoftPion = false;
      double softPiMass = 0.14542; // pion mass in D*->D0pi decay
      double lowMassLimitSoftPion = softPiMass - 3.0 * sigma;
      double highMassLimitSoftPion = softPiMass + 3.0 * sigma;

      if (isD0Cand) {
        if (mDstar1 - part2.mLambda() > 0.) {
          mHistogramRegistry->fill(HIST("SoftPion/softPionMassVsPt"), mDstar1 - part2.mLambda(), part2.pt());
        }
        if ((std::abs(mDstar1 - part2.mLambda()) > lowMassLimitSoftPion) && (std::abs(mDstar1 - part2.mLambda()) < highMassLimitSoftPion)) {
          isSoftPion = true;
        }
      }

      if (isD0barCand) {
        if (mDstar2 - part2.mAntiLambda() > 0.) {
          mHistogramRegistry->fill(HIST("SoftPion/softPionMassVsPt"), mDstar2 - part2.mAntiLambda(), part2.pt());
        }
        if ((std::abs(mDstar2 - part2.mAntiLambda()) > lowMassLimitSoftPion) && (std::abs(mDstar2 - part2.mAntiLambda()) < highMassLimitSoftPion)) {
          isSoftPion = true;
        }
      }
      return isSoftPion;
    } else {
      LOG(fatal) << "FemtoUniverseSoftPionRemoval: Combination of objects not defined - quitting!";
      return false;
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                ///< For QA output
  static constexpr o2::aod::femtouniverseparticle::ParticleType kPartOneType = partOne; ///< Type of particle 1
  static constexpr o2::aod::femtouniverseparticle::ParticleType kPartTwoType = partTwo; ///< Type of particle 2
};
} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESOFTPIONREMOVAL_H_
