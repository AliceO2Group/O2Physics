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

/// \file FastJetUtilities.h
/// \brief Jet related utilities that require fastjet
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_FASTJETUTILITIES_H_
#define PWGJE_CORE_FASTJETUTILITIES_H_

#include <CommonConstants/PhysicsConstants.h>

#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <vector>

constexpr int invalidIndex = -99999999;

enum class JetConstituentStatus {
  invalidStatus = -1,
  track = 0,
  cluster = 1,
  candidate = 2
};

namespace fastjetutilities
{

// Class defined to store additional info which is passed to the FastJet object
class fastjet_user_info : public fastjet::PseudoJet::UserInfoBase
{
  JetConstituentStatus status; // the status of each particle (Options are: TrueParticle (final state particles in generator event which arent special), HFParticle (heavy-flavour particle of interest in generator event), ThermalParticle (particles belonging to the thermal backgound), DecaySisterParticle (other particles poduced in the decay resulting in a non-prompt heavy-flavour particle of interest))
  int index;  // a number unique to each particle in the event

 public:
  fastjet_user_info()
  {
    status = JetConstituentStatus::invalidStatus;
    index = -9;
  }
  fastjet_user_info(JetConstituentStatus _status, int _index)
  {
    status = _status;
    index = _index;
  }
  ~fastjet_user_info() = default;
  void setStatus(JetConstituentStatus set) { status = set; }
  void setIndex(int set) { index = set; }
  JetConstituentStatus getStatus() const { return status; }
  int getIndex() const { return index; }
};

/**
 * Set the fastjet_user_info object when filling the jet constituents.
 *
 * @param constituents vector of constituents to be clustered.
 * @param index global index of constituent
 * @param status status of constituent type
 */

void setFastJetUserInfo(std::vector<fastjet::PseudoJet>& constituents, int index = invalidIndex, JetConstituentStatus status = JetConstituentStatus::track);

/**
 * Add track as a pseudojet object to the fastjet vector
 *
 * @param constituent constituent to be added
 * @param constituents vector of constituents
 * @param index global index of constituent
 * @param status status of constituent type
 * @param status mass hypothesis for constituent
 */

template <typename T>
void fillTracks(const T& constituent, std::vector<fastjet::PseudoJet>& constituents, int index = invalidIndex, JetConstituentStatus status = JetConstituentStatus::track, float mass = o2::constants::physics::MassPiPlus)
{
  if (status == JetConstituentStatus::track || status == JetConstituentStatus::candidate) {
    // auto p = std::sqrt((constituent.px() * constituent.px()) + (constituent.py() * constituent.py()) + (constituent.pz() * constituent.pz()));
    auto energy = std::sqrt((constituent.p() * constituent.p()) + (mass * mass));
    constituents.emplace_back(constituent.px(), constituent.py(), constituent.pz(), energy);
  }
  setFastJetUserInfo(constituents, index, status);
}

/**
 * Add cluster as a pseudojet object to the fastjet vector
 *
 * @param constituent constituent to be added
 * @param constituents vector of constituents
 * @param index global index of constituent
 * @param status status of constituent type
 */

template <typename T>
void fillClusters(const T& constituent, std::vector<fastjet::PseudoJet>& constituents, int index = invalidIndex, int hadronicCorrectionType = 0, JetConstituentStatus status = JetConstituentStatus::cluster)
{
  if (status == JetConstituentStatus::cluster) {
    float constituentEnergy = 0.0;
    if (hadronicCorrectionType == 0) {
      constituentEnergy = constituent.energy();
    }
    if (hadronicCorrectionType == 1) {
      constituentEnergy = constituent.energyCorrectedOneTrack1();
    }
    if (hadronicCorrectionType == 2) {
      constituentEnergy = constituent.energyCorrectedOneTrack2();
    }
    if (hadronicCorrectionType == 3) {
      constituentEnergy = constituent.energyCorrectedAllTracks1();
    }
    if (hadronicCorrectionType == 4) {
      constituentEnergy = constituent.energyCorrectedAllTracks2();
    }
    float constituentPt = constituentEnergy / std::cosh(constituent.eta());
    constituents.emplace_back(constituentPt * std::cos(constituent.phi()), constituentPt * std::sin(constituent.phi()), constituentPt * std::sinh(constituent.eta()), constituentEnergy);
  }
  setFastJetUserInfo(constituents, index, status);
}

}; // namespace fastjetutilities

#endif // PWGJE_CORE_FASTJETUTILITIES_H_
