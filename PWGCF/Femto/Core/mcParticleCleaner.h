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

/// \file tripletCleaner.h
/// \brief triplet cleaner class
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_MCPARTICLECLEANER_H_
#define PWGCF_FEMTO_CORE_MCPARTICLECLEANER_H_

#include "Framework/Configurable.h"

#include <vector>

namespace o2::analysis::femto
{
namespace mcparticlecleaner
{

template <const char* Prefix>
struct ConfMcParticleCleaner : o2::framework::ConfigurableGroup {
  std::string prefix = std::string(Prefix);
  o2::framework::Configurable<bool> cutOnMcTruth{"cutOnMcTruth", false, "Perform cuts based on MC truth information. If false, all following options are ignored"};
  o2::framework::Configurable<int> pdgCode{"pdgCode", 221, "Only consider particles with this exact pdg code (including the sign!)"};
  o2::framework::Configurable<bool> rejectedParticleWithoutMcInformation{"rejectedParticleWithoutMcInformation", true, "If true, all particles which have no associated MC information, are rejected by default"};
  o2::framework::Configurable<std::vector<int>> requiredMotherPdgCodes{"requiredMotherPdgCodes", {}, "Only consider particles whose mothers have one of the supplied pdg codes (inclduing the sign!)"};
  o2::framework::Configurable<std::vector<int>> rejectMotherPdgCodes{"rejectMotherPdgCodes", {}, "Only consider particles whose mothers do not have one of the supplied pdg codes (inclduing the sign!)"};
  o2::framework::Configurable<std::vector<int>> requiredPartonicMotherPdgCodes{"requiredPartonicMotherPdgCodes", {}, "Only consider particles whose partonic mothers have one of the supplied pdg codes (inclduing the sign!)"};
  o2::framework::Configurable<std::vector<int>> rejectPartonicMotherPdgCodes{"rejectPartonicMotherPdgCodes", {}, "Only consider particles whose mothers do not have one of the supplied pdg codes (inclduing the sign!)"};
};

class McParticleCleaner
{
 public:
  McParticleCleaner() = default;
  ~McParticleCleaner() = default;

  template <typename T>
  void init(T const& confMpc)
  {
    mCutOnMcTruth = confMpc.cutOnMcTruth.value;
    mRejectParticleWithoutMcInformation = confMpc.rejectedParticleWithoutMcInformation.value;
    mPdgCode = confMpc.pdgCode.value;
    mRequiredMotherPdgCodes = confMpc.requiredMotherPdgCodes.value;
    mRejectMotherPdgCodes = confMpc.rejectMotherPdgCodes.value;
    mRequiredPartonicMotherPdgCodes = confMpc.requiredPartonicMotherPdgCodes.value;
    mRejectPartonicMotherPdgCodes = confMpc.rejectPartonicMotherPdgCodes.value;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool checkCuts(T1 const& particle, T2 const& /*mcParticles*/, T3 const& /*mcMothers*/, T4 const& /*mcPartonicMothers*/)
  {
    // check whether we apply cuts at all
    if (!mCutOnMcTruth) {
      return true;
    }
    // check if we even have an associated mc particle
    if (!particle.has_fMcParticle()) {
      if (mRejectParticleWithoutMcInformation) {
        return false;
      } else {
        return true;
      }
    }
    // perfrom cuts based on mc information of the particle itself
    auto mcParticle = particle.template fMcParticle_as<T2>();
    if (mPdgCode != mcParticle.pdgCode()) {
      return false;
    }

    // perfrom cuts based on mc information of the mothers
    auto mother = particle.template fMcMother_as<T3>();

    // if list is empty, set it to true and skip the looop
    bool hasMotherWithRequiredPdgCode = mRequiredMotherPdgCodes.empty();
    for (int pdgCode : mRequiredMotherPdgCodes) {
      if (pdgCode == mother.pdgCode()) {
        hasMotherWithRequiredPdgCode = true;
        break;
      }
    }

    bool hasMotherWithRejectedPdgCode = true;
    for (int pdgCode : mRejectMotherPdgCodes) {
      if (pdgCode == mother.pdgCode()) {
        hasMotherWithRejectedPdgCode = false;
        break;
      }
    }

    // perfrom cuts based on mc information of the partonic mothers
    auto partonicMother = particle.template fMcPartMoth_as<T4>();

    // if list is empty, set it to true and skip the looop
    bool hasPartonicMotherWithRequiredPdgCode = mRequiredPartonicMotherPdgCodes.empty();
    for (int pdgCode : mRequiredPartonicMotherPdgCodes) {
      if (pdgCode == partonicMother.pdgCode()) {
        hasPartonicMotherWithRequiredPdgCode = true;
        break;
      }
    }

    bool hasPartonicMotherWithRejectedPdgCode = true;
    for (int pdgCode : mRejectPartonicMotherPdgCodes) {
      if (pdgCode == partonicMother.pdgCode()) {
        hasPartonicMotherWithRejectedPdgCode = false;
        break;
      }
    }

    return hasMotherWithRequiredPdgCode && !hasMotherWithRejectedPdgCode &&
           hasPartonicMotherWithRequiredPdgCode && !hasPartonicMotherWithRejectedPdgCode;
  }

 private:
  bool mCutOnMcTruth = false;
  bool mRejectParticleWithoutMcInformation = true;
  int mPdgCode = 0;
  std::vector<int> mRequiredMotherPdgCodes{};
  std::vector<int> mRejectMotherPdgCodes{};
  std::vector<int> mRequiredPartonicMotherPdgCodes{};
  std::vector<int> mRejectPartonicMotherPdgCodes{};
};

} // namespace mcparticlecleaner
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_MCPARTICLECLEANER_H_
