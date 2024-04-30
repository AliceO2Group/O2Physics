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

template <typename TMCParticle, typename TMCParticles>
int IsFromBeauty(TMCParticle const& p, TMCParticles const& mcparticles)
{
  if (!p.has_mothers()) {
    return -999;
  }

  int motherid = p.mothersIds()[0]; // first mother index
  while (motherid > -1) {
    if (motherid < mcparticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcparticles.iteratorAt(motherid);
      if (abs(mp.pdgCode()) < 1e+9 && (std::to_string(abs(mp.pdgCode()))[std::to_string(abs(mp.pdgCode())).length() - 3] == '5' || std::to_string(abs(mp.pdgCode()))[std::to_string(abs(mp.pdgCode())).length() - 4] == '5')) {
        return motherid;
      }
      if (mp.has_mothers()) {
        motherid = mp.mothersIds()[0];
      } else {
        return -999;
      }
    } else {
      LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid, mcparticles.size());
    }
  }

  return -999;
}

template <typename TMCParticle, typename TMCParticles>
int IsFromCharm(TMCParticle const& p, TMCParticles const& mcparticles)
{
  if (!p.has_mothers()) {
    return -999;
  }

  int motherid = p.mothersIds()[0]; // first mother index
  while (motherid > -1) {
    if (motherid < mcparticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcparticles.iteratorAt(motherid);
      if (abs(mp.pdgCode()) < 1e+9 && (std::to_string(abs(mp.pdgCode()))[std::to_string(abs(mp.pdgCode())).length() - 3] == '4' || std::to_string(abs(mp.pdgCode()))[std::to_string(abs(mp.pdgCode())).length() - 4] == '4')) {
        return motherid;
      }
      if (mp.has_mothers()) {
        motherid = mp.mothersIds()[0];
      } else {
        return -999;
      }
    } else {
      LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid, mcparticles.size());
    }
  }

  return -999;
}

template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
int IsHF(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
{
  if (!p1.has_mothers() || !p2.has_mothers()) {
    return static_cast<int>(EM_HFeeType::kUndef);
  }

  if (p1.mothersIds()[0] == p2.mothersIds()[0]) { // reject same mother. e.g. jspi 443
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

  bool is_direct_from_b1 = abs(mothers_pdg1[0]) < 1e+9 && (std::to_string(mothers_pdg1[0])[std::to_string(mothers_pdg1[0]).length() - 3] == '5' || std::to_string(mothers_pdg1[0])[std::to_string(mothers_pdg1[0]).length() - 4] == '5');
  bool is_direct_from_b2 = abs(mothers_pdg2[0]) < 1e+9 && (std::to_string(mothers_pdg2[0])[std::to_string(mothers_pdg2[0]).length() - 3] == '5' || std::to_string(mothers_pdg2[0])[std::to_string(mothers_pdg2[0]).length() - 4] == '5');
  bool is_prompt_c1 = abs(mothers_pdg1[0]) < 1e+9 && (std::to_string(mothers_pdg1[0])[std::to_string(mothers_pdg1[0]).length() - 3] == '4' || std::to_string(mothers_pdg1[0])[std::to_string(mothers_pdg1[0]).length() - 4] == '4') && IsFromBeauty(p1, mcparticles) < 0;
  bool is_prompt_c2 = abs(mothers_pdg2[0]) < 1e+9 && (std::to_string(mothers_pdg2[0])[std::to_string(mothers_pdg2[0]).length() - 3] == '4' || std::to_string(mothers_pdg2[0])[std::to_string(mothers_pdg2[0]).length() - 4] == '4') && IsFromBeauty(p2, mcparticles) < 0;
  bool is_c_from_b1 = abs(mothers_pdg1[0]) < 1e+9 && (std::to_string(mothers_pdg1[0])[std::to_string(mothers_pdg1[0]).length() - 3] == '4' || std::to_string(mothers_pdg1[0])[std::to_string(mothers_pdg1[0]).length() - 4] == '4') && IsFromBeauty(p1, mcparticles) > 0;
  bool is_c_from_b2 = abs(mothers_pdg2[0]) < 1e+9 && (std::to_string(mothers_pdg2[0])[std::to_string(mothers_pdg2[0]).length() - 3] == '4' || std::to_string(mothers_pdg2[0])[std::to_string(mothers_pdg2[0]).length() - 4] == '4') && IsFromBeauty(p2, mcparticles) > 0;

  if (is_prompt_c1 && is_prompt_c2 && p1.pdgCode() * p2.pdgCode() < 0) {
    mothers_id1.clear();
    mothers_pdg1.clear();
    mothers_id2.clear();
    mothers_pdg2.clear();
    mothers_id1.shrink_to_fit();
    mothers_pdg1.shrink_to_fit();
    mothers_id2.shrink_to_fit();
    mothers_pdg2.shrink_to_fit();
    return static_cast<int>(EM_HFeeType::kCe_Ce); // cc->ee, decay type = 0
  } else if (is_direct_from_b1 && is_direct_from_b2 && p1.pdgCode() * p2.pdgCode() < 0) {
    mothers_id1.clear();
    mothers_pdg1.clear();
    mothers_id2.clear();
    mothers_pdg2.clear();
    mothers_id1.shrink_to_fit();
    mothers_pdg1.shrink_to_fit();
    mothers_id2.shrink_to_fit();
    mothers_pdg2.shrink_to_fit();
    return static_cast<int>(EM_HFeeType::kBe_Be); // bb->ee, decay type = 2
  } else if (is_c_from_b1 && is_c_from_b2 && p1.pdgCode() * p2.pdgCode() < 0) {
    mothers_id1.clear();
    mothers_pdg1.clear();
    mothers_id2.clear();
    mothers_pdg2.clear();
    mothers_id1.shrink_to_fit();
    mothers_pdg1.shrink_to_fit();
    mothers_id2.shrink_to_fit();
    mothers_pdg2.shrink_to_fit();
    return static_cast<int>(EM_HFeeType::kBCe_BCe); // b->c->e and b->c->e, decay type = 1
  } else if ((is_direct_from_b1 && is_c_from_b2) || (is_direct_from_b2 && is_c_from_b1)) {
    if (p1.pdgCode() * p2.pdgCode() < 0) { // ULS
      for (auto& mid1 : mothers_id1) {
        for (auto& mid2 : mothers_id2) {
          if (mid1 == mid2) {
            auto common_mp = mcparticles.iteratorAt(mid1);
            int mp_pdg = common_mp.pdgCode();
            if (abs(mp_pdg) < 1e+9 && (std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 3] == '5' || std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 4] == '5')) {
              mothers_id1.clear();
              mothers_pdg1.clear();
              mothers_id2.clear();
              mothers_pdg2.clear();
              mothers_id1.shrink_to_fit();
              mothers_pdg1.shrink_to_fit();
              mothers_id2.shrink_to_fit();
              mothers_pdg2.shrink_to_fit();
              return static_cast<int>(EM_HFeeType::kBCe_Be_SameB); // b->c->e and b->e, decay type = 3. this should happen only in ULS.
            }
          }
        }    // end of motherid2
      }      // end of motherid1
    } else { // LS
      bool is_same_mother_found = false;
      for (auto& mid1 : mothers_id1) {
        for (auto& mid2 : mothers_id2) {
          if (mid1 == mid2) {
            auto common_mp = mcparticles.iteratorAt(mid1);
            int mp_pdg = common_mp.pdgCode();
            if (abs(mp_pdg) < 1e+9 && (std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 3] == '5' || std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 4] == '5')) {
              is_same_mother_found = true;
            }
          }
        } // end of motherid2
      }   // end of motherid1
      if (!is_same_mother_found) {
        mothers_id1.clear();
        mothers_pdg1.clear();
        mothers_id2.clear();
        mothers_pdg2.clear();
        mothers_id1.shrink_to_fit();
        mothers_pdg1.shrink_to_fit();
        mothers_id2.shrink_to_fit();
        mothers_pdg2.shrink_to_fit();
        return static_cast<int>(EM_HFeeType::kBCe_Be_DiffB); // b->c->e and b->e, decay type = 4. this should happen only in LS. But, this may happen, when ele/pos is reconstructed as pos/ele wrongly. and create LS pair
      }
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
