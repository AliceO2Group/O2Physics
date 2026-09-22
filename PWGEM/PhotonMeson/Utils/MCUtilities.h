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

#include "PWGEM/PhotonMeson/Utils/ParticleOrigin.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/Concepts.h>

#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <ranges>
#include <vector>

//_______________________________________________________________________
namespace o2::aod::pwgem::photonmeson::utils::mcutil
{

constexpr float kVertexEps = 1e-4f; // cm

template <o2::soa::is_iterator TTrack>
bool IsPhysicalPrimary(TTrack const& mctrack)
{
  // This is to check mctrack is ALICE physical primary.
  return (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator());
}
//_______________________________________________________________________
template <o2::soa::is_iterator TCollision, o2::soa::is_iterator T, o2::soa::is_table TMCs>
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
template <o2::soa::is_iterator T, o2::soa::is_table TMCs>
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
template <o2::soa::is_iterator T, o2::soa::is_table TMCs, std::ranges::input_range TTargetPDGs>
int FindMotherInChain(T const& mcparticle, TMCs const& mcparticles, TTargetPDGs const& motherpdgs, const int Depth = 50) // o2-linter: disable=pdg/explicit-code (false positive)
{
  if (!mcparticle.has_mothers() || Depth < 1) {
    return -1;
  }

  int motherid = mcparticle.mothersIds()[0];
  auto mother = mcparticles.iteratorAt(motherid);
  if (std::find(motherpdgs.begin(), motherpdgs.end(), mother.pdgCode()) != motherpdgs.end()) {
    return mcparticle.globalIndex(); // The mother has the required pdg code, so return its daughters global mc particle code.
  }
  return FindMotherInChain(mother, mcparticles, motherpdgs, Depth - 1);
}

//_______________________________________________________________________
/// \brief Go up the decay chain of a mcparticle looking for a mother with the given pdg codes,
/// and return that MOTHER's own index (unlike FindMotherInChain, which returns its daughter).
/// Two different particles that share the same meson ancestor will resolve to the same value here.
/// \param mcparticle iterator of McParticles
/// \param mcparticles table of McParticles
/// \param motherpdgs ranges of mother PDG values to compare against
/// \param Depth how many links should this go up
template <o2::soa::is_iterator T, o2::soa::is_table TMCs, std::ranges::input_range TTargetPDGs>
int GetMesonInChain(T const& mcparticle, TMCs const& mcparticles, TTargetPDGs const& motherpdgs, const int Depth = 15)
{
  int decayChildIdx = FindMotherInChain(mcparticle, mcparticles, motherpdgs, Depth);
  if (decayChildIdx < 0) {
    return -1;
  }
  auto decayChild = mcparticles.iteratorAt(decayChildIdx);
  return decayChild.mothersIds()[0]; // the meson itself, not its daughter
}

//_______________________________________________________________________
/// \brief Go up the decay chain of a mcparticle looking for a mother with the given pdg codes, if found return this mothers daughter
/// E.g. Find the gamma that was created in a pi0 or eta decay
/// \param mcIter iterator of mcparticle -- WILL BE MODIFIED/CONSUMED by this function
/// \param motherPdgs target mother PDG values
/// \param depth how many steps in the chain this check should go maximum before failing
template <o2::soa::is_iterator T, std::ranges::input_range TTargetPDGs>
int findMotherInChain(T& mcIter, TTargetPDGs const& motherPdgs, const int depth = 50)
{
  int currentIndex = mcIter.globalIndex(); // the node whose immediate mother we're about to test

  for (int d = 0; d < depth; ++d) {
    if (!mcIter.has_mothers()) {
      return -1;
    }
    const int motherId = mcIter.mothersIds()[0];
    mcIter.setCursor(motherId);
    if (std::find(motherPdgs.begin(), motherPdgs.end(), mcIter.pdgCode()) != motherPdgs.end()) {
      return currentIndex; // mother matches -- return the node directly below it
    }
    currentIndex = motherId; // no match -- this mother becomes "current" for the next step up
  }
  return -1;
}

//_______________________________________________________________________
template <o2::soa::is_iterator T, o2::soa::is_table TMCs>
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
template <o2::soa::is_iterator TMCParticle, o2::soa::is_table TMCParticles, std::ranges::input_range TTargetPDGs>
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
  return is_equal;
}
//_______________________________________________________________________
template <o2::soa::is_iterator TMCParticle, o2::soa::is_table TMCParticles, std::ranges::input_range TTargetPDGs>
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
  return is_equal;
}
//_______________________________________________________________________
template <o2::soa::is_iterator TMCPhoton, o2::soa::is_table TMCParticles>
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
template <o2::soa::is_iterator TMCParticle, o2::soa::is_table TMCParticles>
bool isGammaGammaDecay(TMCParticle const& mcParticle, TMCParticles const& mcParticles)
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
/// \brief Go up the decay chain of a mcparticle looking for a mother with the given pdg codes, if found return true else false
/// E.g. if electron cluster is coming from a photon return true, if primary electron return false
/// \param mcparticle iterator of mxparticle, WILL BE CHANGED by this function!
/// \param motherPDG target mother PDG value
/// \param depth how many steps in the chain this check should go maximum before failing
template <o2::soa::is_iterator T>
bool isMotherPDG(T& mcparticle, const int motherPDG, const int depth = 10) // o2-linter: disable=pdg/explicit-code (false positive)
{
  if (!mcparticle.has_mothers() || depth < 1) {
    return false;
  }

  int motherid = mcparticle.mothersIds()[0];
  mcparticle.setCursor(motherid);
  if (mcparticle.pdgCode() == motherPDG) {
    return true; // The mother has the required pdg code, so return its daughters global mc particle code.
  }
  return isMotherPDG(mcparticle, motherPDG, depth - 1);
}

//_______________________________________________________________________
/// \brief Go up the decay chain of a mcparticle looking for a mother with the given pdg codes, if found return true else false
/// E.g. if electron cluster is coming from a photon return true, if primary electron return false
/// \param mcparticle iterator of mcparticle, NOT modified by this function
/// \param mcparticleWorking a second iterator of the SAME table, used as scratch space to walk up the chain -- caller must supply this so the function doesn't construct its own
/// \param motherPDG target mother PDG value
/// \param depth how many steps in the chain this check should go maximum before failing
template <o2::soa::is_iterator T>
bool isMotherPDG(const T& mcparticle, T& mcparticleWorking, const int motherPDG, const int depth = 10) // o2-linter: disable=pdg/explicit-code (false positive)
{
  if (!mcparticle.has_mothers() || depth < 1) {
    return false;
  }

  int motherid = mcparticle.mothersIds()[0];
  mcparticleWorking.setCursor(motherid);
  if (mcparticleWorking.pdgCode() == motherPDG) {
    return true; // The mother has the required pdg code.
  }
  return isMotherPDG(mcparticleWorking, mcparticleWorking, motherPDG, depth - 1);
}

//_______________________________________________________________________
/// \brief Check if given particle is from Bremsstrahlung
/// \param mcCursor iterator of mcparticle
/// \param iter shared iterator used to walk to the mother
template <o2::soa::is_iterator TIter>
bool isFromBremsstrahlung(TIter const& mcCursor, TIter& iter, int& motherId)
{
  if (!mcCursor.has_mothers()) {
    return false;
  }
  if (mcCursor.pdgCode() != PDG_t::kGamma) {
    return false; // only a photon can itself be a bremsstrahlung emission
  }
  if (mcCursor.mothersIds().size() != 1) {
    return false; // mother can be only a single lepton, otherwise it might be e+e- annihilation or something else
  }
  motherId = mcCursor.mothersIds()[0];
  iter.setCursor(motherId);
  return std::abs(iter.pdgCode()) == PDG_t::kElectron;
}

//_______________________________________________________________________
/// \brief Go up the decay chain of a mcparticle looking for a mother with the given pdg codes, if found return id else -1
/// E.g. if electron cluster is coming from a photon return true, if primary electron return false
/// \param mcParticle iterator of mxparticle, WILL BE CHANGED by this function!
/// \param motherPDG target mother PDG value
/// \param depth how many steps in the chain this check should go maximum before failing
template <o2::soa::is_iterator T>
int32_t getMotherIndexFromChain(T& mcParticle, const int motherPDG, const int depth = 10)
{
  for (int d = 0; d < depth; ++d) {
    if (!mcParticle.has_mothers()) {
      return -1;
    }
    const int32_t motherid = mcParticle.mothersIds()[0];
    mcParticle.setCursor(motherid);
    if (mcParticle.pdgCode() == motherPDG) {
      return motherid;
    }
  }
  return -1;
}

//_______________________________________________________________________
/// \brief Obtains given photon mcpartlices origin type
/// \param mcPhoton mcparticle iterator of photon
/// \param iter mcparticle iterator used to walk to the mother
/// \param motherDPGs list of pdg values of mothers that would be excepted for PhotonOrigin::Decay
/// \return given photon mcpartlices origin type
template <o2::soa::is_iterator TIter, std::ranges::input_range TTargetPDGs>
o2::analysis::em::PhotonOrigin getPhotonOriginType(TIter const& mcPhoton, TIter& mcIter, TTargetPDGs const& motherPdgs, int& motherId)
{
  switch (mcPhoton.getProcess()) {
    case TMCProcess::kPBrem:
      return o2::analysis::em::PhotonOrigin::Bremsstrahlung;
    case TMCProcess::kPAnnihilation:
      return o2::analysis::em::PhotonOrigin::Annihilation;
    case TMCProcess::kPHadronic:
      return o2::analysis::em::PhotonOrigin::Hadronic;
    default:
      break;
  }
  if (!mcPhoton.has_mothers()) {
    return o2::analysis::em::PhotonOrigin::Other;
  }
  mcIter.setCursor(mcPhoton.globalIndex());
  motherId = findMotherInChain(mcIter, motherPdgs);
  if (motherId >= 0) {
    return o2::analysis::em::PhotonOrigin::Decay;
  }
  motherId = mcPhoton.mothersIds()[0];
  mcIter.setCursor(motherId);
  if ((std::abs(mcIter.pdgCode()) >= PDG_t::kDown && std::abs(mcIter.pdgCode()) <= PDG_t::kTop) || std::abs(mcIter.pdgCode()) == PDG_t::kGluon) {
    return o2::analysis::em::PhotonOrigin::Direct;
  }
  return o2::analysis::em::PhotonOrigin::Other;
}

//_______________________________________________________________________
/// \brief Obtains given lepton mcpartlices origin type
/// \param mcLepton mcparticle iterator of lepton
/// \param iter mcparticle iterator used to walk to the mother
/// \param motherDPGs list of pdg values of mothers that would be excepted for PhotonOrigin::Decay
/// \return given lepton mcpartlices origin type
template <o2::soa::is_iterator TIter, std::ranges::input_range TTargetPDGs>
o2::analysis::em::LeptonOrigin getLeptonOriginType(TIter const& mcLepton, TIter& mcIter, TTargetPDGs const& motherPdgs)
{
  switch (mcLepton.getProcess()) {
    case TMCProcess::kPPair:
      return o2::analysis::em::LeptonOrigin::Conversion;
    case TMCProcess::kPCompton:
      return o2::analysis::em::LeptonOrigin::Compton;
    case TMCProcess::kPPhotoelectric:
      return o2::analysis::em::LeptonOrigin::PhotoElectric;
    case TMCProcess::kPDeltaRay:
      return o2::analysis::em::LeptonOrigin::DeltaRay;
    case TMCProcess::kPDecay: {
      if (!mcLepton.has_mothers()) {
        return o2::analysis::em::LeptonOrigin::Other;
      }
      mcIter.setCursor(mcLepton.mothersIds()[0]);
      if (std::find(motherPdgs.begin(), motherPdgs.end(), mcIter.pdgCode()) != motherPdgs.end()) {
        return o2::analysis::em::LeptonOrigin::DirectMesonDecay;
      }
      return o2::analysis::em::LeptonOrigin::Other; // decay, but not from your target meson list
    }
    default:
      return o2::analysis::em::LeptonOrigin::Other;
  }
}

//_______________________________________________________________________
/// \brief Photon mother of a conversion electron, addressed by MC index
/// \param mcParticles table of MC particles (aod::McParticles or aod::EMMCParticles)
/// \param mcId index of the MC particle behind a track or V0 leg, -1 if the track carries no label
/// \return index of the mother photon if the particle is a transport electron with a photon as first mother, else -1
template <o2::soa::is_table TMCParticles>
int photonMotherId(TMCParticles const& mcParticles, const int mcId)
{
  if (mcId < 0) {
    return -1;
  }
  return IsEleFromPC(mcParticles.iteratorAt(mcId), mcParticles);
}

//_______________________________________________________________________
/// \brief MC leg information for a V0 candidate
struct PhotonMCInfo {
  int mcPosId = -1;               /// MC particle of the positive leg, -1 = no label
  int mcNegId = -1;               /// MC particle of the negative leg, -1 = no label
  int posPhotonId = -1;           /// photon mother of the positive leg (photonMotherId), -1 = not a conversion electron
  int negPhotonId = -1;           /// same for the negative leg
  int motherId = -1;              /// common first mother of both legs, whatever its species, -1 = none
  int motherPdg = 0;              /// PDG code of that common mother
  bool isTruePhoton = false;      /// both legs are conversion electrons of the SAME photon
  bool isPhysicalPrimary = false; /// the common mother is a physical primary
  [[nodiscard]] bool hasBothLabels() const { return mcPosId >= 0 && mcNegId >= 0; }
  [[nodiscard]] bool bothLegsFromPhotons() const { return posPhotonId >= 0 && negPhotonId >= 0; }
};

//_______________________________________________________________________
/// \brief Fills PhotonMCInfo from the two leg labels
/// \param mcParticles table of MC particles
/// \param mcPosId MC index of the positive leg, -1 if unlabelled
/// \param mcNegId MC index of the negative leg, -1 if unlabelled
/// \return truth record of the candidate; with one unlabelled leg only that leg's fields stay at -1
template <o2::soa::is_table TMCParticles>
PhotonMCInfo makePhotonMCInfo(TMCParticles const& mcParticles, const int mcPosId, const int mcNegId)
{
  PhotonMCInfo info;
  info.mcPosId = mcPosId;
  info.mcNegId = mcNegId;
  info.posPhotonId = photonMotherId(mcParticles, mcPosId);
  info.negPhotonId = photonMotherId(mcParticles, mcNegId);
  info.isTruePhoton = (info.posPhotonId >= 0) && (info.posPhotonId == info.negPhotonId);
  if (!info.hasBothLabels()) {
    return info;
  }
  const auto mcPos = mcParticles.iteratorAt(mcPosId);
  const auto mcNeg = mcParticles.iteratorAt(mcNegId);
  const int motherPos = mcPos.has_mothers() ? mcPos.mothersIds()[0] : -1;
  const int motherNeg = mcNeg.has_mothers() ? mcNeg.mothersIds()[0] : -1;
  if (motherPos < 0 || motherPos != motherNeg) {
    return info; // no common valid mother
  }
  const auto mother = mcParticles.iteratorAt(motherPos);
  info.motherId = motherPos;
  info.motherPdg = mother.pdgCode();
  info.isPhysicalPrimary = mother.isPhysicalPrimary();
  return info;
}

//_______________________________________________________________________
/// \brief Fills PhotonMCInfo for a V0 photon of the skimmed tables
/// \tparam TLegs leg table joined with the MC labels (e.g. soa::Join<aod::V0Legs, aod::V0LegMCLabels>)
/// \param v0 iterator of the V0 photon
/// \param mcParticles table of aod::EMMCParticles
/// \return truth record of this V0
template <typename TLegs, o2::soa::is_iterator TV0, o2::soa::is_table TMCParticles>
PhotonMCInfo makePhotonMCInfo(TV0 const& v0, TMCParticles const& mcParticles)
{
  const auto pos = v0.template posTrack_as<TLegs>();
  const auto neg = v0.template negTrack_as<TLegs>();
  return makePhotonMCInfo(mcParticles,
                          pos.has_emmcparticle() ? static_cast<int>(pos.emmcparticleId()) : -1,
                          neg.has_emmcparticle() ? static_cast<int>(neg.emmcparticleId()) : -1);
}

//_______________________________________________________________________
/// \brief Truth class of a pair of two V0 candidates
enum class PairTruthType : int {
  Unknown = 0,
  TrueTrueDistinct,   /// two different true photons
  TrueTrueSamePhoton, /// the same true photon twice (duplicate V0)
  SharedMcLeg,        /// the two V0s contain the same MC particle as a leg
  TrueFake,           /// one true photon, one fake
  FakeFake,           /// two fakes
  Pi0Daughters,       /// two true photons from the same pi0
  NClasses
};

//_______________________________________________________________________
/// \brief Truth class of a pair from the two V0 records, without the pi0 check
/// \param m1 truth record of the first V0
/// \param m2 truth record of the second V0
/// \return PairTruthType; TrueTrueSamePhoton has priority over SharedMcLeg, the shared-leg test runs across charges
[[nodiscard]] inline PairTruthType classifyPairTruth(PhotonMCInfo const& m1, PhotonMCInfo const& m2)
{
  if (m1.isTruePhoton && m2.isTruePhoton && m1.posPhotonId == m2.posPhotonId) {
    return PairTruthType::TrueTrueSamePhoton;
  }
  const bool sharedLeg =
    (m1.mcPosId >= 0 && (m1.mcPosId == m2.mcPosId || m1.mcPosId == m2.mcNegId)) ||
    (m1.mcNegId >= 0 && (m1.mcNegId == m2.mcPosId || m1.mcNegId == m2.mcNegId));
  if (sharedLeg) {
    return PairTruthType::SharedMcLeg;
  }
  if (!m1.isTruePhoton && !m2.isTruePhoton) {
    return PairTruthType::FakeFake;
  }
  if (m1.isTruePhoton != m2.isTruePhoton) {
    return PairTruthType::TrueFake;
  }
  return PairTruthType::TrueTrueDistinct;
}

//_______________________________________________________________________
/// \brief Are two true photons daughters of the same pi0
/// \param m1 truth record of the first V0
/// \param m2 truth record of the second V0
/// \param mcParticles table of MC particles
/// \return true if both V0s are true photons and IsXFromY finds the same pi0 as first mother of both photons
template <o2::soa::is_table TMCParticles>
bool isPi0DaughterPair(PhotonMCInfo const& m1, PhotonMCInfo const& m2, TMCParticles const& mcParticles)
{
  if (!m1.isTruePhoton || !m2.isTruePhoton) {
    return false;
  }
  const int pi0Of1 = IsXFromY(mcParticles.iteratorAt(m1.posPhotonId), mcParticles, PDG_t::kGamma, PDG_t::kPi0);
  const int pi0Of2 = IsXFromY(mcParticles.iteratorAt(m2.posPhotonId), mcParticles, PDG_t::kGamma, PDG_t::kPi0);
  return pi0Of1 >= 0 && pi0Of1 == pi0Of2;
}

//_______________________________________________________________________
/// \brief Full truth class of a pair, as the HBT task bins it
/// \param m1 truth record of the first V0
/// \param m2 truth record of the second V0
/// \param mcParticles table of MC particles
/// \return classifyPairTruth, with TrueTrueDistinct promoted to Pi0Daughters where isPi0DaughterPair holds
template <o2::soa::is_table TMCParticles>
PairTruthType pairTruthType(PhotonMCInfo const& m1, PhotonMCInfo const& m2, TMCParticles const& mcParticles)
{
  const PairTruthType t = classifyPairTruth(m1, m2);
  if (t == PairTruthType::TrueTrueDistinct && isPi0DaughterPair(m1, m2, mcParticles)) {
    return PairTruthType::Pi0Daughters;
  }
  return t;
}

//_______________________________________________________________________
/// \brief Classifies the V0 cross-wise built candidates. The crossed pairs are evaluted whether it is from the same photon or not.
/// \param m1 truth record of the first V0
/// \param m2 truth record of the second V0
/// \return 1 if both directions cross (full 2x2 cross: same two photons, legs swapped), 2 if only one direction crosses, 0 otherwise
[[nodiscard]] inline int classifyFakeSubtype(PhotonMCInfo const& m1, PhotonMCInfo const& m2)
{
  const bool cross12 = m1.posPhotonId >= 0 && m1.posPhotonId == m2.negPhotonId;
  const bool cross21 = m2.posPhotonId >= 0 && m2.posPhotonId == m1.negPhotonId;
  if (cross12 && cross21) {
    return 1;
  }
  return (cross12 || cross21) ? 2 : 0; // o2-linter: disable=magic-number (class code, see brief)
}

//_______________________________________________________________________
/// \brief Particle a mother chain ends at
enum class AncestorKind : int {
  NoneFound = 0,   /// no common ancestor within the search depth
  StringOrCluster, /// PYTHIA string or cluster, or a bare parton (quark, gluon)
  BeamOrNucleus,   /// beam particle, nucleus, nucleon or diquark remnant
  Pi0,
  Eta,
  OtherMeson,
  Baryon,
  Photon, // the ancestor is a photon: one particle descends from the other
  Electron,
  Other,
  NClasses
};

//_______________________________________________________________________
/// \brief Maps a PDG code onto AncestorKind
/// \param pdgIn PDG code, the sign is ignored
/// \return AncestorKind of the code; ranges of codes, not single particles
inline AncestorKind ancestorKindOfPdg(int pdgIn)
{
  const int pdg = std::abs(pdgIn);
  constexpr int kStringLo = 91, kStringHi = 94;
  constexpr int kDiquarkLo = 1103, kDiquarkHi = 5503;
  constexpr int kNucleusOffset = 1000000000;
  constexpr int kMesonHi = 999;
  constexpr int kBaryonHi = 9999;
  constexpr int kQuarkHi = 8, kGluon = 21;
  if ((pdg >= kStringLo && pdg <= kStringHi) || (pdg >= 1 && pdg <= kQuarkHi) || pdg == kGluon) {
    return AncestorKind::StringOrCluster;
  }
  if (pdg >= kNucleusOffset || (pdg >= kDiquarkLo && pdg <= kDiquarkHi)) {
    return AncestorKind::BeamOrNucleus;
  }
  if (pdg == PDG_t::kPi0) {
    return AncestorKind::Pi0;
  }
  if (pdg == o2::constants::physics::Pdg::kEta) {
    return AncestorKind::Eta;
  }
  if (pdg == PDG_t::kGamma) {
    return AncestorKind::Photon;
  }
  if (pdg == PDG_t::kElectron) {
    return AncestorKind::Electron;
  }
  if (pdg == PDG_t::kProton || pdg == PDG_t::kNeutron) {
    return AncestorKind::BeamOrNucleus;
  }
  if (pdg > 0 && pdg <= kMesonHi) {
    return AncestorKind::OtherMeson;
  }
  if (pdg > kMesonHi && pdg <= kBaryonHi) {
    return AncestorKind::Baryon;
  }
  return AncestorKind::Other;
}

//_______________________________________________________________________
/// \brief The particle a chain originates from, stopping before the beam
/// \param mcParticles table of MC particles
/// \param id index of the particle to start from
/// \param maxGen safety cap against self-referencing mother links, not a physics limit
/// \return index of the origin particle, or -1 if id was invalid
template <o2::soa::is_table TMCParticles>
int originParticleId(TMCParticles const& mcParticles, int id, int maxGen = 50)
{
  if (id < 0) {
    return -1;
  }
  int cur = id;
  for (int gen = 0; gen < maxGen; ++gen) {
    const auto p = mcParticles.iteratorAt(cur);
    if (!p.has_mothers()) {
      return cur;
    }
    const int mother = p.mothersIds()[0];
    if (mother < 0 || mother == cur) {
      return cur;
    }
    const AncestorKind kindOfMother = ancestorKindOfPdg(mcParticles.iteratorAt(mother).pdgCode());
    if (kindOfMother == AncestorKind::BeamOrNucleus || kindOfMother == AncestorKind::StringOrCluster) {
      return cur;
    }
    cur = mother;
  }
  return cur;
}

//_______________________________________________________________________
/// \brief Chains meetup like in findCommonAncestor
struct CommonAncestor {
  int id{-1};   /// index of the lowest common ancestor, -1 if none was found
  int gen1{-1}; /// generations from the first particle up to it (0 = it IS that particle)
  int gen2{-1}; /// the same for the second particle
};

//_______________________________________________________________________
/// \brief Lowest common ancestor of two MC particles
/// \param mcParticles table of MC particles
/// \param id1 index of the first particle
/// \param id2 index of the second particle
/// \param maxGen how many generations to walk on each side
/// \param stopAtBeam treat beam, nucleus and string as "no mother"
/// \return CommonAncestor with id = -1 if the chains do not meet within maxGen
template <o2::soa::is_table TMCParticles>
CommonAncestor findCommonAncestor(TMCParticles const& mcParticles, int id1, int id2,
                                  int maxGen = 4, bool stopAtBeam = false)
{
  CommonAncestor out;
  if (id1 < 0 || id2 < 0) {
    return out;
  }
  auto motherOf = [&mcParticles, stopAtBeam](int id) -> int {
    const auto particle = mcParticles.iteratorAt(id);
    if (!particle.has_mothers()) {
      return -1;
    }
    const int mother = particle.mothersIds()[0];
    if (mother < 0 || mother == id) {
      return -1;
    }
    if (stopAtBeam) {
      const AncestorKind kind = ancestorKindOfPdg(mcParticles.iteratorAt(mother).pdgCode());
      if (kind == AncestorKind::BeamOrNucleus || kind == AncestorKind::StringOrCluster) {
        return -1;
      }
    }
    return mother;
  };
  std::vector<int> chain1;
  chain1.reserve(static_cast<size_t>(maxGen) + 1);
  for (int id = id1, gen = 0; id >= 0 && gen <= maxGen; ++gen) {
    chain1.push_back(id);
    id = motherOf(id);
  }
  for (int id = id2, gen = 0; id >= 0 && gen <= maxGen; ++gen) {
    const auto it = std::find(chain1.begin(), chain1.end(), id);
    if (it != chain1.end()) {
      out.id = id;
      out.gen1 = static_cast<int>(std::distance(chain1.begin(), it));
      out.gen2 = gen;
      return out;
    }
    id = motherOf(id);
  }
  return out;
}

//_______________________________________________________________________
/// \brief Class of the lowest common ancestor of two MC particles
/// \param mcParticles table of MC particles
/// \param id1 index of the first particle
/// \param id2 index of the second particle
/// \param maxGen how many generations to walk on each side
/// \return AncestorKind of the common ancestor, AncestorKind::NoneFound when the chains do not meet within maxGen
template <o2::soa::is_table TMCParticles>
AncestorKind commonAncestorKind(TMCParticles const& mcParticles, int id1, int id2, int maxGen = 4)
{
  const CommonAncestor anc = findCommonAncestor(mcParticles, id1, id2, maxGen, false);
  if (anc.id < 0) {
    return AncestorKind::NoneFound;
  }
  return ancestorKindOfPdg(mcParticles.iteratorAt(anc.id).pdgCode());
}

//_______________________________________________________________________
/// \brief Distinct parents for mc truth
///   V0_a = (e+_1, e-_1), V0_b = (e+_2, e-_2)   N=2, consistent=2, shared=no   signal
///   V0_a = (e+_1, e-_2), V0_b = (e+_2, e-_1)   N=2, consistent=0, shared=yes  leg swap
///   V0_a = (e+_1, e-_1), V0_b = (e+_2, e-_3)   N=3, consistent=1, shared=no   one real + one swap
///   V0_a = (e+_1, e-_2), V0_b = (e+_3, e-_4)   N=4, consistent=0, shared=no   pure combinatorics
///   V0_a = (e+_1, e-_1), V0_b = (e+_1',e-_1')  N=1, consistent=2, shared=yes  duplicate (split track)
///   V0_a = (e+_1, X),    V0_b = (Y, e-_1)      N=1, consistent=0, shared=yes  one photon, two foreign tracks
struct MotherCensus {
  int nConvLegs = 0;   /// legs whose mother is a photon (0-4)
  int nMothers = 0;    /// distinct parent photons among them (0-4)
  int nConsistent = 0; /// V0s whose two legs share the same parent photon (0-2)
  bool shared = false; /// one parent photon appears in BOTH V0s
  std::array<int, 4> distinct{{-1, -1, -1, -1}};
  [[nodiscard]] int splitCode() const { return 2 * nConsistent + (shared ? 1 : 0); } // o2-linter: disable=magic-number (two values per nConsistent)
};

//_______________________________________________________________________
/// \brief Census of the four legs of a photon pair
/// \param motherPhotonIds the parent photon of each leg, in the order (V0 A e+, V0 A e-, V0 B e+, V0 B e-);
/// \return MotherCensus of the four legs
inline MotherCensus censusMothers(std::array<int, 4> const& motherPhotonIds)
{
  constexpr int kNLegs = 4;
  const auto& mp = motherPhotonIds;
  MotherCensus c;
  for (int i = 0; i < kNLegs; ++i) {
    const auto ui = static_cast<size_t>(i);
    if (mp[ui] < 0) {
      continue;
    }
    c.nConvLegs++;
    bool seen = false;
    for (int j = 0; j < i; ++j) { // only the first appearance counts as a new mother
      seen = seen || (mp[static_cast<size_t>(j)] == mp[ui]);
    }
    if (!seen) {
      c.distinct[static_cast<size_t>(c.nMothers)] = mp[ui]; // nMothers < 4: at most four legs
      c.nMothers++;
    }
  }
  c.nConsistent = ((mp[0] >= 0 && mp[0] == mp[1]) ? 1 : 0) + ((mp[2] >= 0 && mp[2] == mp[3]) ? 1 : 0);
  for (int i = 0; i < 2; ++i) { // o2-linter: disable=magic-number (the two legs of the first V0)
    for (int j = 2; j < kNLegs; ++j) {
      c.shared = c.shared || (mp[static_cast<size_t>(i)] >= 0 && mp[static_cast<size_t>(i)] == mp[static_cast<size_t>(j)]);
    }
  }
  return c;
}

//_______________________________________________________________________
/// \brief Amount of V0 mothers per legs.
/// \param m1 truth record of the first V0
/// \param m2 truth record of the second V0
/// \return censusMothers of (m1 e+, m1 e-, m2 e+, m2 e-)
[[nodiscard]] inline MotherCensus censusMothers(PhotonMCInfo const& m1, PhotonMCInfo const& m2)
{
  return censusMothers({m1.posPhotonId, m1.negPhotonId, m2.posPhotonId, m2.negPhotonId});
}

//_______________________________________________________________________
/// \brief Analysis for two photons based on their MC track label content.
enum class DupClass : int {
  Unknown = 0,  /// at least one candidate without both labels
  SameMcPhoton, /// a duplicate
  SharedMcLeg,  /// share an MC particle as a leg, but different mother photons
  CrossBuilt,   /// legs of two true photons built crosswise (swap)
  DistinctTrue, /// two different, cleanly reconstructed photons
  InvolvesFake, /// at least one is not a true photon, no swap
  NClasses
};

//_______________________________________________________________________
/// \brief Duplicate class of two candidates from their truth records
/// \param a truth record of the first candidate
/// \param b truth record of the second candidate
/// \return DupClass; the shared-leg test runs before the label requirement, a V0 with one labelled leg can still share it
[[nodiscard]] inline DupClass classifyDupPair(PhotonMCInfo const& a, PhotonMCInfo const& b)
{
  const bool sharedLeg =
    (a.mcPosId >= 0 && (a.mcPosId == b.mcPosId || a.mcPosId == b.mcNegId)) ||
    (a.mcNegId >= 0 && (a.mcNegId == b.mcPosId || a.mcNegId == b.mcNegId));
  if (!a.hasBothLabels() || !b.hasBothLabels()) {
    return sharedLeg ? DupClass::SharedMcLeg : DupClass::Unknown;
  }
  if (a.isTruePhoton && b.isTruePhoton && a.posPhotonId == b.posPhotonId) {
    return DupClass::SameMcPhoton;
  }
  if (sharedLeg) {
    return DupClass::SharedMcLeg;
  }
  if (classifyFakeSubtype(a, b) != 0) {
    return DupClass::CrossBuilt;
  }
  if (a.isTruePhoton && b.isTruePhoton) {
    return DupClass::DistinctTrue;
  }
  return DupClass::InvolvesFake;
}

//_______________________________________________________________________
/// \brief How one V0 candidate was built, by MC truth
enum class V0RecoStatus : int {
  NoLabel = 0, ///< at least one leg without MC label
  NoConvLeg,   ///< no leg is a conversion electron
  OneConvLeg,  ///< one conversion electron, one foreign track
  CrossBuilt,  ///< two conversion electrons of DIFFERENT photons
  TruePhoton,  ///< two conversion electrons of the same photon
  NClasses
};

//_______________________________________________________________________
/// \brief Reco status of a candidate from its truth record
/// \param m truth record of the candidate
/// \return V0RecoStatus
[[nodiscard]] inline V0RecoStatus classifyV0Reco(PhotonMCInfo const& m)
{
  if (!m.hasBothLabels()) {
    return V0RecoStatus::NoLabel;
  }
  if (m.isTruePhoton) {
    return V0RecoStatus::TruePhoton;
  }
  if (m.bothLegsFromPhotons()) {
    return V0RecoStatus::CrossBuilt;
  }
  return (m.posPhotonId >= 0 || m.negPhotonId >= 0) ? V0RecoStatus::OneConvLeg : V0RecoStatus::NoConvLeg;
}

//_______________________________________________________________________
/// \brief True conversion point of a photon: the production vertex of its first e+/e- daughter
/// \param mcParticles table of MC particles
/// \param photonId index of the photon
/// \param xyz output, (vx, vy, vz) of the conversion (cm)
/// \return true if the photon has an electron daughter, false otherwise (xyz untouched)
template <o2::soa::is_table TMCParticles>
bool conversionPoint(TMCParticles const& mcParticles, const int photonId, std::array<float, 3>& xyz)
{
  if (photonId < 0) {
    return false;
  }
  const auto photon = mcParticles.iteratorAt(photonId);
  if (!photon.has_daughters()) {
    return false;
  }
  for (const auto& dId : photon.daughtersIds()) {
    if (dId < 0) {
      continue;
    }
    const auto d = mcParticles.iteratorAt(dId);
    if (std::abs(d.pdgCode()) == PDG_t::kElectron) {
      xyz = {static_cast<float>(d.vx()), static_cast<float>(d.vy()), static_cast<float>(d.vz())};
      return true;
    }
  }
  return false;
}

//_______________________________________________________________________
/// \brief True q_inv of two MC photons (massless)
/// \param mcParticles table of MC particles
/// \param id1 index of the first photon
/// \param id2 index of the second photon
/// \return q_inv = sqrt(2 (E1 E2 - p1.p2)) in GeV/c, -1 if an index is invalid
template <o2::soa::is_table TMCParticles>
float trueQinvOfPhotons(TMCParticles const& mcParticles, const int id1, const int id2)
{
  if (id1 < 0 || id2 < 0) {
    return -1.f;
  }
  const auto pa = mcParticles.iteratorAt(id1);
  const auto pb = mcParticles.iteratorAt(id2);
  const float e1 = std::hypot(pa.px(), pa.py(), pa.pz());
  const float e2 = std::hypot(pb.px(), pb.py(), pb.pz());
  return std::sqrt(std::max(0.f, 2.f * (e1 * e2 - pa.px() * pb.px() - pa.py() * pb.py() - pa.pz() * pb.pz())));
}

//_______________________________________________________________________
/// \brief True opening angle of two MC photons
/// \param mcParticles table of MC particles
/// \param id1 index of the first photon
/// \param id2 index of the second photon
/// \return angle between the two momenta in rad, -1 if an index is invalid
template <o2::soa::is_table TMCParticles>
float trueOpeningAngle(TMCParticles const& mcParticles, const int id1, const int id2)
{
  if (id1 < 0 || id2 < 0) {
    return -1.f;
  }
  const auto pa = mcParticles.iteratorAt(id1);
  const auto pb = mcParticles.iteratorAt(id2);
  const float n1 = std::hypot(pa.px(), pa.py(), pa.pz());
  const float n2 = std::hypot(pb.px(), pb.py(), pb.pz());
  if (n1 < 1e-9f || n2 < 1e-9f) { // o2-linter: disable=magic-number (zero momentum guard)
    return -1.f;
  }
  const float c = (pa.px() * pb.px() + pa.py() * pb.py() + pa.pz() * pb.pz()) / (n1 * n2);
  return std::acos(std::clamp(c, -1.f, 1.f));
}

} // namespace o2::aod::pwgem::photonmeson::utils::mcutil

#endif // PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_
