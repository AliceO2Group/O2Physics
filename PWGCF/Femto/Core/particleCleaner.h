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

/// \file particleCleaner.h
/// \brief particle cleaner class
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PARTICLECLEANER_H_
#define PWGCF_FEMTO_CORE_PARTICLECLEANER_H_

#include "Framework/Configurable.h"

#include <string>
#include <vector>

namespace o2::analysis::femto
{
namespace particlecleaner
{

template <const char* Prefix>
struct ConfParticleCleaner : o2::framework::ConfigurableGroup {
  std::string prefix = std::string(Prefix);
  o2::framework::Configurable<bool> activate{"activate", false, "Activate particle cleaner"};
  o2::framework::Configurable<std::vector<int>> requiredPdgCodes{"requiredPdgCodes", {}, "Only consider particles with this exact pdg code (including the sign!)"};
  o2::framework::Configurable<std::vector<int>> rejectedPdgCodes{"rejectedPdgCodes", {}, "Reject particles with this exact pdg code (including the sign!)"};
  o2::framework::Configurable<bool> rejectParticleWithoutMcParticle{"rejectParticleWithoutMcParticle", false, "If true, particles which have no associated MC information, are rejected"};
  o2::framework::Configurable<bool> rejectParticleWithoutMcMother{"rejectParticleWithoutMcMother", false, "If true, particles which have no associated mother are rejected"};
  o2::framework::Configurable<bool> rejectParticleWithoutMcPartonicMother{"rejectParticleWithoutMcPartonicMother", false, "If true, all particles which have no associated partonic mother"};
  o2::framework::Configurable<std::vector<int>> requiredMotherPdgCodes{"requiredMotherPdgCodes", {}, "Only consider particles whose mothers have one of the supplied pdg codes (inclduing the sign!)"};
  o2::framework::Configurable<std::vector<int>> rejectMotherPdgCodes{"rejectMotherPdgCodes", {}, "Only consider particles whose mothers do not have one of the supplied pdg codes (inclduing the sign!)"};
  o2::framework::Configurable<std::vector<int>> requiredPartonicMotherPdgCodes{"requiredPartonicMotherPdgCodes", {}, "Only consider particles whose partonic mothers have one of the supplied pdg codes (inclduing the sign!)"};
  o2::framework::Configurable<std::vector<int>> rejectPartonicMotherPdgCodes{"rejectPartonicMotherPdgCodes", {}, "Only consider particles whose mothers do not have one of the supplied pdg codes (inclduing the sign!)"};
};

constexpr const char PrefixTrackCleaner1[] = "TrackCleaner1";
constexpr const char PrefixTrackCleaner2[] = "TrackCleaner2";
constexpr const char PrefixTrackCleaner3[] = "TrackCleaner3";
using ConfTrackCleaner1 = ConfParticleCleaner<PrefixTrackCleaner1>;
using ConfTrackCleaner2 = ConfParticleCleaner<PrefixTrackCleaner2>;
using ConfTrackCleaner3 = ConfParticleCleaner<PrefixTrackCleaner3>;

constexpr const char PrefixLambdaCleaner1[] = "LambdaCleaner1";
constexpr const char PrefixLambdaCleaner2[] = "LambdaCleaner2";
using ConfLambdaCleaner1 = ConfParticleCleaner<PrefixLambdaCleaner1>;
using ConfLambdaCleaner2 = ConfParticleCleaner<PrefixLambdaCleaner2>;

constexpr const char PrefixK0shortCleaner1[] = "K0shortCleaner1";
constexpr const char PrefixK0shortCleaner2[] = "K0shortCleaner2";
using ConfK0shortCleaner1 = ConfParticleCleaner<PrefixK0shortCleaner1>;
using ConfK0shortCleaner2 = ConfParticleCleaner<PrefixK0shortCleaner2>;

constexpr const char PrefixSigmaCleaner1[] = "SigmaCleaner1";
constexpr const char PrefixSigmaCleaner2[] = "SigmaCleaner2";
using ConfSigmaCleaner1 = ConfParticleCleaner<PrefixSigmaCleaner1>;
using ConfSigmaCleaner2 = ConfParticleCleaner<PrefixSigmaCleaner2>;

constexpr const char PrefixSigmaPlusCleaner1[] = "SigmaPlusCleaner1";
constexpr const char PrefixSigmaPlusCleaner2[] = "SigmaPlusCleaner2";
using ConfSigmaPlusCleaner1 = ConfParticleCleaner<PrefixSigmaPlusCleaner1>;
using ConfSigmaPlusCleaner2 = ConfParticleCleaner<PrefixSigmaPlusCleaner2>;

constexpr const char PrefixXiCleaner1[] = "XiCleaner1";
constexpr const char PrefixXiCleaner2[] = "XiCleaner2";
using ConfXiCleaner1 = ConfParticleCleaner<PrefixXiCleaner1>;
using ConfXiCleaner2 = ConfParticleCleaner<PrefixXiCleaner2>;

constexpr const char PrefixOmegaCleaner1[] = "OmegaCleaner1";
constexpr const char PrefixOmegaCleaner2[] = "OmegaCleaner2";
using ConfOmegaCleaner1 = ConfParticleCleaner<PrefixOmegaCleaner1>;
using ConfOmegaCleaner2 = ConfParticleCleaner<PrefixOmegaCleaner2>;

class ParticleCleaner
{
 public:
  ParticleCleaner() = default;
  ~ParticleCleaner() = default;

  template <typename T>
  void init(T const& confMpc)
  {
    mActivate = confMpc.activate.value;

    mRejectParticleWithoutMcParticle = confMpc.rejectParticleWithoutMcParticle.value;
    mRejectParticleWithoutMcMother = confMpc.rejectParticleWithoutMcMother.value;
    mRejectParticleWithoutMcPartonicMother = confMpc.rejectParticleWithoutMcPartonicMother.value;

    mRequiredPdgCodes = confMpc.requiredPdgCodes.value;
    mRejectedPdgCodes = confMpc.rejectedPdgCodes.value;

    mRequiredMotherPdgCodes = confMpc.requiredMotherPdgCodes.value;
    mRejectedMotherPdgCodes = confMpc.rejectMotherPdgCodes.value;

    mRequiredPartonicMotherPdgCodes = confMpc.requiredPartonicMotherPdgCodes.value;
    mRejectedPartonicMotherPdgCodes = confMpc.rejectPartonicMotherPdgCodes.value;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool isClean(T1 const& particle,
               T2 const& /*mcParticles*/,
               T3 const& /*mcMothers*/,
               T4 const& /*mcPartonicMothers*/)
  {
    if (!mActivate) {
      return true;
    }

    bool hasRequiredPdgCode = true;
    bool hasRejectedPdgCode = false;

    bool hasMotherWithRequiredPdgCode = true;
    bool hasMotherWithRejectedPdgCode = false;

    bool hasPartonicMotherWithRequiredPdgCode = true;
    bool hasPartonicMotherWithRejectedPdgCode = false;

    // MC particle
    if (!particle.has_fMcParticle()) {
      if (mRejectParticleWithoutMcParticle || !mRequiredPdgCodes.empty()) {
        return false;
      }
    } else {
      auto mcParticle = particle.template fMcParticle_as<T2>();

      if (!mRequiredPdgCodes.empty()) {
        hasRequiredPdgCode = false;
        for (int const& pdgCode : mRequiredPdgCodes) {
          if (pdgCode == mcParticle.pdgCode()) {
            hasRequiredPdgCode = true;
            break;
          }
        }
      }

      if (!mRejectedPdgCodes.empty()) {
        for (int const& pdgCode : mRejectedPdgCodes) {
          if (pdgCode == mcParticle.pdgCode()) {
            hasRejectedPdgCode = true;
            break;
          }
        }
      }
    }

    // MC mother
    if (!particle.has_fMcMother()) {
      if (mRejectParticleWithoutMcMother || !mRequiredMotherPdgCodes.empty()) {
        return false;
      }
    } else {
      auto mother = particle.template fMcMother_as<T3>();

      if (!mRequiredMotherPdgCodes.empty()) {
        hasMotherWithRequiredPdgCode = false;
        for (int const& pdgCode : mRequiredMotherPdgCodes) {
          if (pdgCode == mother.pdgCode()) {
            hasMotherWithRequiredPdgCode = true;
            break;
          }
        }
      }

      if (!mRejectedMotherPdgCodes.empty()) {
        for (int const& pdgCode : mRejectedMotherPdgCodes) {
          if (pdgCode == mother.pdgCode()) {
            hasMotherWithRejectedPdgCode = true;
            break;
          }
        }
      }
    }

    // MC partonic mother
    if (!particle.has_fMcPartMoth()) {
      if (mRejectParticleWithoutMcPartonicMother ||
          !mRequiredPartonicMotherPdgCodes.empty()) {
        return false;
      }
    } else {
      auto partonicMother = particle.template fMcPartMoth_as<T4>();

      if (!mRequiredPartonicMotherPdgCodes.empty()) {
        hasPartonicMotherWithRequiredPdgCode = false;
        for (int const& pdgCode : mRequiredPartonicMotherPdgCodes) {
          if (pdgCode == partonicMother.pdgCode()) {
            hasPartonicMotherWithRequiredPdgCode = true;
            break;
          }
        }
      }

      if (!mRejectedPartonicMotherPdgCodes.empty()) {
        for (int const& pdgCode : mRejectedPartonicMotherPdgCodes) {
          if (pdgCode == partonicMother.pdgCode()) {
            hasPartonicMotherWithRejectedPdgCode = true;
            break;
          }
        }
      }
    }

    return hasRequiredPdgCode && !hasRejectedPdgCode &&
           hasMotherWithRequiredPdgCode && !hasMotherWithRejectedPdgCode &&
           hasPartonicMotherWithRequiredPdgCode &&
           !hasPartonicMotherWithRejectedPdgCode;
  }

 private:
  bool mActivate = false;
  bool mRejectParticleWithoutMcParticle = true;
  bool mRejectParticleWithoutMcMother = true;
  bool mRejectParticleWithoutMcPartonicMother = true;
  std::vector<int> mRequiredPdgCodes = {};
  std::vector<int> mRejectedPdgCodes = {};
  std::vector<int> mRequiredMotherPdgCodes{};
  std::vector<int> mRejectedMotherPdgCodes{};
  std::vector<int> mRequiredPartonicMotherPdgCodes{};
  std::vector<int> mRejectedPartonicMotherPdgCodes{};
};

} // namespace particlecleaner
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_PARTICLECLEANER_H_
