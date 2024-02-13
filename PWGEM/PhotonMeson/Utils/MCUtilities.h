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

/// \commonly used for MC analysis.
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_

#include <vector>
#include <algorithm>
#include "Framework/AnalysisTask.h"

//_______________________________________________________________________
namespace o2::aod::pwgem::mcutil
{
template <typename TCollision, typename TTrack, typename TMCs>
bool IsPhysicalPrimary(TCollision const& mccollision, TTrack const& mctrack, TMCs const& mcTracks)
{
  // This is to check mctrack is ALICE physical primary.
  // https://inspirehep.net/files/4c26ef5fb432df99bdc1ff847653502f

  if (mctrack.isPhysicalPrimary()) { // this is the first priority. In fact, this does not happen to neutral mesons in ALICE.
    return true;
  }

  if (!mctrack.producedByGenerator()) {
    return false;
  }
  float r3D = sqrt(pow(mctrack.vx() - mccollision.posX(), 2) + pow(mctrack.vy() - mccollision.posY(), 2) + pow(mctrack.vz() - mccollision.posZ(), 2)); // cm
  if (r3D > 1.0) {
    return false;
  }

  // exclude weak decay. K0S and Lambda are the 2 most relevant strange particles decaying into neutral mesons.
  if (mctrack.has_mothers()) {
    // auto mp = mctrack.template mothers_first_as<TMCs>();
    int motherid = mctrack.mothersIds()[0]; // first mother index
    while (motherid > -1) {
      if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
        auto mp = mcTracks.iteratorAt(motherid);
        int pdg_mother = mp.pdgCode();
        // LOGF(info, "mctrack.globalIndex() = %d, mp.globalIndex() = %d , pdg_mother = %d", mctrack.globalIndex(), mp.globalIndex(), pdg_mother);
        if (abs(pdg_mother) == 310 || abs(pdg_mother) == 130 || abs(pdg_mother) == 3122) {
          return false;
        }
        if (mp.has_mothers()) {
          motherid = mp.mothersIds()[0]; // first mother index
          auto mp_tmp = mcTracks.iteratorAt(motherid);
          int pdg_mother_tmp = mp_tmp.pdgCode();
          if (pdg_mother_tmp == pdg_mother) { // strange protection. mother of 111 is 111. observed in LHC23k6d on 25.January.2024
            return false;
          }
        } else {
          motherid = -999;
        }
      }
    }

    // for (auto& m : mctrack.mothersIds()) {
    //   if (m < mcTracks.size()) { // protect against bad mother indices
    //     auto mp = mcTracks.iteratorAt(m);
    //     int pdg_mother = mp.pdgCode();
    //     if (abs(pdg_mother) == 310 || abs(pdg_mother) == 3122 || abs(pdg_mother) == 3212) {
    //       return false;
    //     }
    //   }
    // }
  }
  return true;
}
//_______________________________________________________________________
template <typename TCollision, typename T, typename TMCs>
bool IsFromWD(TCollision const& mccollision, T const& mctrack, TMCs const& mcTracks)
{
  // is this particle from weak decay?

  if (mctrack.isPhysicalPrimary()) { // this is the first priority.
    return false;
  }

  if (mctrack.producedByGenerator()) {
    return false;
  }

  if (mctrack.has_mothers()) {
    // auto mp = mctrack.template mothers_first_as<TMCs>();
    int motherid = mctrack.mothersIds()[0]; // first mother index
    while (motherid > -1) {
      if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
        auto mp = mcTracks.iteratorAt(motherid);
        int pdg_mother = mp.pdgCode();
        if (abs(pdg_mother) == 310 || abs(pdg_mother) == 130 || abs(pdg_mother) == 3122) {
          // LOGF(info, "mctrack.globalIndex() = %d, mp.globalIndex() = %d , pdg_mother = %d", mctrack.globalIndex(), mp.globalIndex(), pdg_mother);
          return true;
        }
        if (mp.has_mothers()) {
          motherid = mp.mothersIds()[0]; // first mother index
        } else {
          motherid = -999;
        }
      }
    }
  } else {
    return false;
  }
  return false;
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
template <typename T, typename TMCs>
int IsEleFromPC(T const& mctrack, TMCs const& mcTracks)
{
  // is election from photon conversion? returns index of mother photon
  if (abs(mctrack.pdgCode()) != 11)
    return -1;
  if (mctrack.producedByGenerator())
    return -1;
  if (mctrack.has_mothers()) {
    int motherid = mctrack.mothersIds()[0]; // first mother
    auto mp = mcTracks.iteratorAt(motherid);
    int pdg_mother = mp.pdgCode();
    if (pdg_mother == 22) {
      return motherid;
    }
  } else {
    return -1;
  }
  return -1;
}
//_______________________________________________________________________
template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
int FindCommonMotherFrom2Prongs(TMCParticle1 const& p1, TMCParticle2 const& p2, const int expected_pdg1, const int expected_pdg2, const int expected_mother_pdg, TMCParticles const& mcparticles)
{
  if (p1.globalIndex() == p2.globalIndex())
    return -1; // mc particle p1 and p2 is identical. reject.

  if (p1.pdgCode() != expected_pdg1)
    return -1;
  if (p2.pdgCode() != expected_pdg2)
    return -1;

  if (!p1.has_mothers())
    return -1;
  if (!p2.has_mothers())
    return -1;

  // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

  int motherid1 = p1.mothersIds()[0];
  auto mother1 = mcparticles.iteratorAt(motherid1);
  int mother1_pdg = mother1.pdgCode();

  int motherid2 = p2.mothersIds()[0];
  auto mother2 = mcparticles.iteratorAt(motherid2);
  int mother2_pdg = mother2.pdgCode();

  // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

  if (motherid1 != motherid2)
    return -1;
  if (mother1_pdg != mother2_pdg)
    return -1;
  if (mother1_pdg != expected_mother_pdg)
    return -1;
  return motherid1;
}
template <typename TMCParticle1, typename TMCParticle2, typename TMCParticle3, typename TMCParticles>
int FindCommonMotherFrom3Prongs(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticle3 const& p3, const int expected_pdg1, const int expected_pdg2, const int expected_pdg3, const int expected_mother_pdg, TMCParticles const& mcparticles)
{
  if (p1.globalIndex() == p2.globalIndex())
    return -1; // mc particle p1 and p2 are identical. reject.
  if (p2.globalIndex() == p3.globalIndex())
    return -1; // mc particle p2 and p3 are identical. reject.
  if (p3.globalIndex() == p1.globalIndex())
    return -1; // mc particle p3 and p1 are identical. reject.

  if (p1.pdgCode() != expected_pdg1)
    return -1;
  if (p2.pdgCode() != expected_pdg2)
    return -1;
  if (p3.pdgCode() != expected_pdg3)
    return -1;

  if (!p1.has_mothers())
    return -1;
  if (!p2.has_mothers())
    return -1;
  if (!p3.has_mothers())
    return -1;

  // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

  int motherid1 = p1.mothersIds()[0];
  auto mother1 = mcparticles.iteratorAt(motherid1);
  int mother1_pdg = mother1.pdgCode();

  int motherid2 = p2.mothersIds()[0];
  auto mother2 = mcparticles.iteratorAt(motherid2);
  int mother2_pdg = mother2.pdgCode();

  int motherid3 = p3.mothersIds()[0];
  auto mother3 = mcparticles.iteratorAt(motherid3);
  int mother3_pdg = mother3.pdgCode();

  // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

  if (motherid1 != motherid2)
    return -1;
  if (motherid2 != motherid3)
    return -1;
  if (motherid3 != motherid1)
    return -1;

  if (mother1_pdg != mother2_pdg)
    return -1;
  if (mother2_pdg != mother3_pdg)
    return -1;
  if (mother3_pdg != mother1_pdg)
    return -1;

  if (mother1_pdg != expected_mother_pdg)
    return -1;
  return motherid1;
}
//_______________________________________________________________________
template <typename TMCParticle, typename TMCParticles, typename TTargetPDGs>
bool IsInAcceptance(TMCParticle const& mcparticle, TMCParticles const& mcparticles, TTargetPDGs const& target_pdgs, const float ymin, const float ymax, const float phimin, const float phimax)
{
  if (mcparticle.y() < ymin || ymax < mcparticle.y()) {
    return false; // mother rapidity is out of acceptance
  }
  if (mcparticle.phi() < phimin || phimax < mcparticle.phi()) {
    return false; // mother rapidity is out of acceptance
  }
  auto daughtersIds = mcparticle.daughtersIds(); // always size = 2. first and last index. one should run loop from the first index to the last index.
  if (daughtersIds[0] < 0 || daughtersIds[1] < 0) {
    return false;
  }

  if (daughtersIds[1] - daughtersIds[0] + 1 != static_cast<int>(target_pdgs.size())) {
    return false;
  }
  std::vector<int> pdgs;
  pdgs.reserve(target_pdgs.size());
  for (int idau = daughtersIds[0]; idau <= daughtersIds[1]; idau++) {
    auto daughter = mcparticles.iteratorAt(idau);
    pdgs.emplace_back(daughter.pdgCode());

    if (daughter.eta() < ymin || ymax < daughter.eta()) {
      return false;
    }
    if (daughter.phi() < phimin || phimax < daughter.phi()) {
      return false;
    }
  } // end of daughter loop

  sort(pdgs.begin(), pdgs.end());
  bool is_equal = std::equal(pdgs.cbegin(), pdgs.cend(), target_pdgs.cbegin());
  pdgs.clear();
  pdgs.shrink_to_fit();
  if (!is_equal) {
    return false; // garantee daughter is in acceptance.
  }
  return true;
}
} // namespace o2::aod::pwgem::mcutil
//_______________________________________________________________________
//_______________________________________________________________________
#endif // PWGEM_PHOTONMESON_UTILS_MCUTILITIES_H_
