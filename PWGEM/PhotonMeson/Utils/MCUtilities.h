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

/// \file MCUtilities.h
/// \brief commonly used for MC analysis.
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_

#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

//_______________________________________________________________________
namespace o2::aod::pwgem::photonmeson::utils::mcutil
{
template <typename TTrack>
bool IsPhysicalPrimary(TTrack const& mctrack)
{
  // This is to check mctrack is ALICE physical primary.
  if (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename TCollision, typename T, typename TMCs>
int IsFromWD(TCollision const&, T const& mctrack, TMCs const& mcTracks)
{
  // is this particle from weak decay?
  if (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) {
    return -1;
  }

  if (mctrack.has_mothers()) {
    // auto mp = mctrack.template mothers_first_as<TMCs>();
    int motherid = mctrack.mothersIds()[0]; // first mother index
    while (motherid > -1) {
      if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
        auto mp = mcTracks.iteratorAt(motherid);
        int pdg_mother = mp.pdgCode();
        if (std::abs(pdg_mother) == kK0Short || std::abs(pdg_mother) == kK0Long || std::abs(pdg_mother) == kLambda0) {
          // LOGF(info, "mctrack.globalIndex() = %d, mp.globalIndex() = %d , pdg_mother = %d", mctrack.globalIndex(), mp.globalIndex(), pdg_mother);
          return motherid;
        }
        if (mp.has_mothers()) {
          motherid = mp.mothersIds()[0]; // first mother index
        } else {
          motherid = -999;
        }
      }
    }
  } else {
    return -1;
  }
  return -1;
}
//_______________________________________________________________________
template <typename T, typename TMCs>
int IsXFromY(T const& mctrack, TMCs const& mcTracks, const int pdgX, const int pdgY)
{
  // is photon from pi0? returns index of mother photon
  if (mctrack.pdgCode() != pdgX) {
    return -1;
  }
  if (mctrack.has_mothers()) {
    int motherid = mctrack.mothersIds()[0]; // first mother
    auto mp = mcTracks.iteratorAt(motherid);
    int pdg_mother = mp.pdgCode();
    if (pdg_mother == pdgY) {
      return motherid;
    }
  } else {
    return -1;
  }
  return -1;
}
//_______________________________________________________________________
// Go up the decay chain of a mcparticle looking for a mother with the given pdg codes, if found return this mothers daughter
// E.g. Find the gamma that was created in a pi0 or eta decay
template <typename T, typename TMCs, typename TTargetPDGs>
int FindMotherInChain(T const& mcparticle, TMCs const& mcparticles, TTargetPDGs const& motherpdgs, const int Depth = 50) // o2-linter: disable=pdg/explicit-code (false positive)
{
  if (!mcparticle.has_mothers() || Depth < 1) {
    return -1;
  }

  int motherid = mcparticle.mothersIds()[0];
  auto mother = mcparticles.iteratorAt(motherid);
  if (std::find(motherpdgs.begin(), motherpdgs.end(), mother.pdgCode()) != motherpdgs.end()) {
    return mcparticle.globalIndex(); // The mother has the required pdg code, so return its daughters global mc particle code.
  } else {
    return FindMotherInChain(mother, mcparticles, motherpdgs, Depth - 1);
  }
}
//_______________________________________________________________________
template <typename T, typename TMCs>
int IsEleFromPC(T const& mctrack, TMCs const& mcTracks)
{
  // is election from photon conversion? returns index of mother photon
  if (std::abs(mctrack.pdgCode()) != kElectron) {
    return -1;
  }
  if (mctrack.producedByGenerator()) {
    return -1;
  }
  if (mctrack.has_mothers()) {
    int motherid = mctrack.mothersIds()[0]; // first mother
    auto mp = mcTracks.iteratorAt(motherid);
    int pdg_mother = mp.pdgCode();
    if (pdg_mother == kGamma) {
      return motherid;
    }
  } else {
    return -1;
  }
  return -1;
}
//_______________________________________________________________________
template <typename TMCParticle, typename TMCParticles, typename TTargetPDGs>
bool IsInAcceptanceNonDerived(TMCParticle const& mcparticle, TMCParticles const& mcparticles, TTargetPDGs target_pdgs, const float ymin, const float ymax, const float phimin, const float phimax)
{
  // contents in vector of daughter ID is different.

  if (mcparticle.y() < ymin || ymax < mcparticle.y()) {
    return false; // mother rapidity is out of acceptance
  }
  if (mcparticle.phi() < phimin || phimax < mcparticle.phi()) {
    return false; // mother rapidity is out of acceptance
  }
  // auto daughtersIds = mcparticle.daughtersIds(); // always size = 2. first and last index. one should run loop from the first index to the last index.
  int ndau = mcparticle.daughtersIds()[1] - mcparticle.daughtersIds()[0] + 1;

  if (ndau != static_cast<int>(target_pdgs.size())) {
    return false;
  }
  std::vector<int> pdgs;
  pdgs.reserve(target_pdgs.size());
  for (int daughterId = mcparticle.daughtersIds()[0]; daughterId <= mcparticle.daughtersIds()[1]; ++daughterId) {
    if (daughterId < 0) {
      pdgs.clear();
      pdgs.shrink_to_fit();
      return false;
    }
    auto daughter = mcparticles.iteratorAt(daughterId);
    pdgs.emplace_back(daughter.pdgCode());

    if (daughter.eta() < ymin || ymax < daughter.eta()) {
      pdgs.clear();
      pdgs.shrink_to_fit();
      return false;
    }
    if (daughter.phi() < phimin || phimax < daughter.phi()) {
      pdgs.clear();
      pdgs.shrink_to_fit();
      return false;
    }
  } // end of daughter loop

  sort(target_pdgs.begin(), target_pdgs.end());
  sort(pdgs.begin(), pdgs.end());
  bool is_equal = std::equal(pdgs.cbegin(), pdgs.cend(), target_pdgs.cbegin());
  pdgs.clear();
  pdgs.shrink_to_fit();
  if (!is_equal) {
    return false; // garantee daughter is in acceptance.
  }
  return true;
}
//_______________________________________________________________________
template <typename TMCParticle, typename TMCParticles, typename TTargetPDGs>
bool IsInAcceptance(TMCParticle const& mcparticle, TMCParticles const& mcparticles, TTargetPDGs target_pdgs, const float ymin, const float ymax, const float phimin, const float phimax)
{
  if (mcparticle.y() < ymin || ymax < mcparticle.y()) {
    return false; // mother rapidity is out of acceptance
  }
  if (mcparticle.phi() < phimin || phimax < mcparticle.phi()) {
    return false; // mother rapidity is out of acceptance
  }
  auto daughtersIds = mcparticle.daughtersIds();

  if (daughtersIds.size() != target_pdgs.size()) {
    return false;
  }
  std::vector<int> pdgs;
  pdgs.reserve(target_pdgs.size());
  for (const auto& daughterId : daughtersIds) {
    if (daughterId < 0) {
      pdgs.clear();
      pdgs.shrink_to_fit();
      return false;
    }
    auto daughter = mcparticles.iteratorAt(daughterId);
    pdgs.emplace_back(daughter.pdgCode());

    if (daughter.eta() < ymin || ymax < daughter.eta()) {
      pdgs.clear();
      pdgs.shrink_to_fit();
      return false;
    }
    if (daughter.phi() < phimin || phimax < daughter.phi()) {
      pdgs.clear();
      pdgs.shrink_to_fit();
      return false;
    }
  } // end of daughter loop

  sort(target_pdgs.begin(), target_pdgs.end());
  sort(pdgs.begin(), pdgs.end());
  bool is_equal = std::equal(pdgs.cbegin(), pdgs.cend(), target_pdgs.cbegin());
  pdgs.clear();
  pdgs.shrink_to_fit();
  if (!is_equal) {
    return false; // garantee daughter is in acceptance.
  }
  return true;
}
//_______________________________________________________________________
template <typename TMCPhoton, typename TMCParticles>
bool IsConversionPointInAcceptance(TMCPhoton const& mcphoton, const float max_r_gen, const float max_eta_gen, const float margin_z_mc, TMCParticles const& mcparticles)
{
  if (std::abs(mcphoton.pdgCode()) != kGamma) {
    return false;
  }

  auto daughtersIds = mcphoton.daughtersIds();
  if (daughtersIds.size() != 2) { // o2-linter: disable=magic-number (2 is not that magic in this context)
    return false;
  }

  for (const auto& daughterId : daughtersIds) {
    if (daughterId < 0) {
      return false;
    }
    auto daughter = mcparticles.iteratorAt(daughterId);
    if (std::abs(daughter.pdgCode()) != kElectron) {
      return false;
    }

    if (daughter.producedByGenerator()) {
      return false;
    }

    float rxy_gen_e = std::sqrt(std::pow(daughter.vx(), 2) + std::pow(daughter.vy(), 2));
    // LOGF(info, "daughterId = %d , pdg = %d , vx = %f , vy = %f , vz = %f, rxy = %f", daughterId, daughter.pdgCode(), daughter.vx(), daughter.vy(), daughter.vz(), rxy_gen_e);
    if (rxy_gen_e > max_r_gen || rxy_gen_e < std::abs(daughter.vz()) * std::tan(2 * std::atan(std::exp(-max_eta_gen))) - margin_z_mc) {
      return false;
    }
  } // end of daughter loop

  return true;
}
//_______________________________________________________________________
template <typename TMCParticle, typename TMCParticles>
bool isGammaGammaDecay(TMCParticle mcParticle, TMCParticles mcParticles)
{
  auto daughtersIds = mcParticle.daughtersIds();
  if (daughtersIds.size() != 2) { // o2-linter: disable=magic-number (2 is not that magic in this context)
    return false;
  }
  for (const auto& daughterId : daughtersIds) {
    if (mcParticles.iteratorAt(daughterId).pdgCode() != kGamma) {
      return false;
    }
  }
  return true;
}
//_______________________________________________________________________
} // namespace o2::aod::pwgem::photonmeson::utils::mcutil
//_______________________________________________________________________
//_______________________________________________________________________
#endif // PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_
