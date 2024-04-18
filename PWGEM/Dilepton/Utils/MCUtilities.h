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

#ifndef PWGEM_DILEPTON_UTILS_MCUTILITIES_H_
#define PWGEM_DILEPTON_UTILS_MCUTILITIES_H_

#include <string>
#include <vector>
#include <algorithm>

//_______________________________________________________________________
namespace o2::aod::pwgem::dilepton::mcutil
{
enum class EM_HFeeType : int {
  kUndef = -1,
  kCe_Ce = 0,        // ULS
  kBe_Be = 1,        // ULS
  kBCe_BCe = 2,      // ULS
  kBCe_Be_SameB = 3, // ULS
  kBCe_Be_DiffB = 4, // LS
};
template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
int IsHF(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
{
  if (!p1.has_mothers() || !p2.has_mothers()) {
    return static_cast<int>(EM_HFeeType::kUndef);
  }

  if (p1.mothersIds()[0] == p2.mothersIds()[0]) { // same mother
    return static_cast<int>(EM_HFeeType::kUndef); // this never happens in correlated HF->ee decays
  }

  // store all mother1 relation
  std::vector<int> mothers_id1;
  std::vector<int> mothers_pdg1;
  int motherid1 = p1.mothersIds()[0]; // first mother index
  while (motherid1 > -1) {
    if (motherid1 < mcparticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcparticles.iteratorAt(motherid1);
      mothers_id1.emplace_back(motherid1);
      mothers_pdg1.emplace_back(mp.pdgCode());

      if (mp.has_mothers()) {
        motherid1 = mp.mothersIds()[0];
      } else {
        motherid1 = -999;
      }
    } else {
      LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid1, mcparticles.size());
    }
  }

  // store all mother2 relation
  std::vector<int> mothers_id2;
  std::vector<int> mothers_pdg2;
  int motherid2 = p2.mothersIds()[0]; // first mother index
  while (motherid2 > -1) {
    if (motherid2 < mcparticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcparticles.iteratorAt(motherid2);
      mothers_id2.emplace_back(motherid2);
      mothers_pdg2.emplace_back(mp.pdgCode());

      if (mp.has_mothers()) {
        motherid2 = mp.mothersIds()[0];
      } else {
        motherid2 = -999;
      }
    } else {
      LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid2, mcparticles.size());
    }
  }

  if (std::to_string(mothers_pdg1[0]).find("5") != std::string::npos && std::to_string(mothers_pdg2[0]).find("5") != std::string::npos) {
    return static_cast<int>(EM_HFeeType::kBe_Be); // bb->ee, decay type = 2
    // this is easy. first mother is b hadron for both leg.
  }

  if (std::to_string(mothers_pdg1[0]).find("4") != std::string::npos && std::to_string(mothers_pdg2[0]).find("4") != std::string::npos) {
    // mother is c hadron. next, check c is prompt or non-prompt.

    bool is_c_from_b1 = false;
    for (unsigned int i1 = 1; i1 < mothers_pdg1.size(); i1++) {
      if (std::to_string(mothers_pdg1[i1]).find("5") != std::string::npos) {
        is_c_from_b1 = true;
        break;
      }
    }
    bool is_c_from_b2 = false;
    for (unsigned int i2 = 1; i2 < mothers_pdg2.size(); i2++) {
      if (std::to_string(mothers_pdg2[i2]).find("5") != std::string::npos) {
        is_c_from_b2 = true;
        break;
      }
    }

    if (!is_c_from_b1 && !is_c_from_b2) {
      return static_cast<int>(EM_HFeeType::kCe_Ce); // prompt cc->ee, decay type = 0
    } else if (is_c_from_b1 && is_c_from_b2) {
      return static_cast<int>(EM_HFeeType::kBCe_BCe); // b->c->e and b->c->e, decay type = 1
    } else {
      for (auto& mid1 : mothers_id1) {
        for (auto& mid2 : mothers_id2) {
          if (mid1 == mid2) {
            return static_cast<int>(EM_HFeeType::kBCe_Be_SameB); // b->c->e and c->e, decay type = 3. this should happen only in ULS.
          }
        }                                                  // end of mother id 2
      }                                                    // end of mother id 1
      return static_cast<int>(EM_HFeeType::kBCe_Be_DiffB); // b->c->e and c->e, decay type = 4. this should happen only in LS. But, this may happen, when ele/pos is reconstructed as pos/ele wrongly. and create LS pair
    }
  }
  mothers_id1.clear();
  mothers_pdg1.clear();
  mothers_id2.clear();
  mothers_pdg2.clear();
  mothers_id1.shrink_to_fit();
  mothers_pdg1.shrink_to_fit();
  mothers_id2.shrink_to_fit();
  mothers_pdg2.shrink_to_fit();
  return static_cast<int>(EM_HFeeType::kUndef);
}
//_______________________________________________________________________
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::mcutil
#endif // PWGEM_DILEPTON_UTILS_MCUTILITIES_H_
