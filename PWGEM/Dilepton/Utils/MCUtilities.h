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
namespace o2::aod::pwgem::dilepton::utils::mcutil
{
enum class EM_HFeeType : int {
  kUndef = -1,
  kCe_Ce = 0,        // ULS
  kBe_Be = 1,        // ULS
  kBCe_BCe = 2,      // ULS
  kBCe_Be_SameB = 3, // ULS
  kBCe_Be_DiffB = 4, // LS
};

//_______________________________________________________________________
template <typename TMCParticle1, typename TMCParticle2>
int FindCommonMotherFrom2ProngsWithoutPDG(TMCParticle1 const& p1, TMCParticle2 const& p2)
{
  if (p1.globalIndex() == p2.globalIndex())
    return -1; // mc particle p1 and p2 is identical. reject.

  if (!p1.has_mothers())
    return -1;
  if (!p2.has_mothers())
    return -1;

  // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

  int motherid1 = p1.mothersIds()[0];
  int motherid2 = p2.mothersIds()[0];

  // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

  if (motherid1 != motherid2)
    return -1;
  return motherid1;
}
//_______________________________________________________________________
template <typename TMCParticle1, typename TMCParticle2, typename TMCParticle3>
int FindCommonMotherFrom3ProngsWithoutPDG(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticle3 const& p3)
{
  if (p1.globalIndex() == p2.globalIndex())
    return -1; // mc particle p1 and p2 are identical. reject.
  if (p2.globalIndex() == p3.globalIndex())
    return -1; // mc particle p2 and p3 are identical. reject.
  if (p3.globalIndex() == p1.globalIndex())
    return -1; // mc particle p3 and p1 are identical. reject.

  if (!p1.has_mothers())
    return -1;
  if (!p2.has_mothers())
    return -1;
  if (!p3.has_mothers())
    return -1;

  // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

  int motherid1 = p1.mothersIds()[0];
  int motherid2 = p2.mothersIds()[0];
  int motherid3 = p3.mothersIds()[0];

  // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

  if (motherid1 != motherid2)
    return -1;
  if (motherid2 != motherid3)
    return -1;
  if (motherid3 != motherid1)
    return -1;
  return motherid1;
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
//_______________________________________________________________________
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
template <typename TMCParticle, typename TMCParticles>
int IsFromBeauty(TMCParticle const& p, TMCParticles const& mcparticles)
{
  if (!p.has_mothers()) {
    return -999;
  }

  int motherid = p.mothersIds()[0]; // first mother index
  auto mp_tmp = mcparticles.iteratorAt(motherid);
  if (abs(mp_tmp.pdgCode()) < 1e+9 && (std::to_string(abs(mp_tmp.pdgCode()))[std::to_string(abs(mp_tmp.pdgCode())).length() - 2] == '5' && std::to_string(abs(mp_tmp.pdgCode()))[std::to_string(abs(mp_tmp.pdgCode())).length() - 3] == '5') && abs(mp_tmp.pdgCode()) % 2 == 1) {
    return -999; // reject bottomonia
  }

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
  auto mp_tmp = mcparticles.iteratorAt(motherid);
  if (abs(mp_tmp.pdgCode()) < 1e+9 && (std::to_string(abs(mp_tmp.pdgCode()))[std::to_string(abs(mp_tmp.pdgCode())).length() - 2] == '4' && std::to_string(abs(mp_tmp.pdgCode()))[std::to_string(abs(mp_tmp.pdgCode())).length() - 3] == '4') && abs(mp_tmp.pdgCode()) % 2 == 1) {
    return -999; // reject bottomonia
  }
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

  bool is_direct_from_b1 = IsFromBeauty(p1, mcparticles) > 0 && IsFromCharm(p1, mcparticles) < 0;
  bool is_direct_from_b2 = IsFromBeauty(p2, mcparticles) > 0 && IsFromCharm(p2, mcparticles) < 0;
  bool is_prompt_c1 = IsFromBeauty(p1, mcparticles) < 0 && IsFromCharm(p1, mcparticles) > 0;
  bool is_prompt_c2 = IsFromBeauty(p2, mcparticles) < 0 && IsFromCharm(p2, mcparticles) > 0;
  bool is_c_from_b1 = IsFromBeauty(p1, mcparticles) > 0 && IsFromCharm(p1, mcparticles) > 0;
  bool is_c_from_b2 = IsFromBeauty(p2, mcparticles) > 0 && IsFromCharm(p2, mcparticles) > 0;

  if (is_direct_from_b1 && is_direct_from_b2 && p1.pdgCode() * p2.pdgCode() < 0) {
    mothers_id1.clear();
    mothers_pdg1.clear();
    mothers_id2.clear();
    mothers_pdg2.clear();
    mothers_id1.shrink_to_fit();
    mothers_pdg1.shrink_to_fit();
    mothers_id2.shrink_to_fit();
    mothers_pdg2.shrink_to_fit();
    return static_cast<int>(EM_HFeeType::kBe_Be); // bb->ee, decay type = 2
  }
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
  }
  if (is_c_from_b1 && is_c_from_b2 && p1.pdgCode() * p2.pdgCode() < 0) {
    mothers_id1.clear();
    mothers_pdg1.clear();
    mothers_id2.clear();
    mothers_pdg2.clear();
    mothers_id1.shrink_to_fit();
    mothers_pdg1.shrink_to_fit();
    mothers_id2.shrink_to_fit();
    mothers_pdg2.shrink_to_fit();
    return static_cast<int>(EM_HFeeType::kBCe_BCe); // b->c->e and b->c->e, decay type = 1
  }
  if ((is_direct_from_b1 && is_c_from_b2) || (is_direct_from_b2 && is_c_from_b1)) {
    if (p1.pdgCode() * p2.pdgCode() < 0) { // ULS
      for (auto& mid1 : mothers_id1) {
        for (auto& mid2 : mothers_id2) {
          if (mid1 == mid2) {
            auto common_mp = mcparticles.iteratorAt(mid1);
            int mp_pdg = common_mp.pdgCode();
            bool is_mp_diquark = (1100 < abs(mp_pdg) && abs(mp_pdg) < 5600) && std::to_string(mp_pdg)[std::to_string(mp_pdg).length() - 2] == '0';
            if (!is_mp_diquark && abs(mp_pdg) < 1e+9 && (std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 3] == '5' || std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 4] == '5')) {
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
            bool is_mp_diquark = (1100 < abs(mp_pdg) && abs(mp_pdg) < 5600) && std::to_string(mp_pdg)[std::to_string(mp_pdg).length() - 2] == '0';
            if (!is_mp_diquark && abs(mp_pdg) < 1e+9 && (std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 3] == '5' || std::to_string(abs(mp_pdg))[std::to_string(abs(mp_pdg)).length() - 4] == '5')) {
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

template <typename T, typename U>
int searchMothers(T& p, U& mcParticles, int pdg, bool equal)
{ // find the first ancestor that is equal/not-equal pdg

  if (!p.has_mothers()) {
    return -1;
  }
  auto mothersids = p.mothersIds();
  std::vector<int> allmothersids;

  if (mothersids.size() == 2) {
    if (mothersids[0] == mothersids[1]) {
      allmothersids.push_back(mothersids[0]);
    } else if (mothersids[1] < mothersids[0]) {
      allmothersids.push_back(mothersids[0]);
      allmothersids.push_back(mothersids[1]);
    } else if ((80 < abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) && abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) < 90) || (100 < abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) && abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) < 110)) { // NOTE: THIS IS GENERATOR DEPENDENT AND WORKS ONLY FOR PYTHIA!
      for (int i = mothersids[0]; i <= mothersids[1]; i++) {
        allmothersids.push_back(i);
      }
    } else {
      allmothersids.push_back(mothersids[0]);
      allmothersids.push_back(mothersids[1]);
    }
  } else {
    allmothersids.push_back(mothersids[0]);
  }

  if (equal) { // we are searching for the quark
    int quark_id = -1;
    int next_mother_id = -1;
    for (int i : allmothersids) {
      auto mother = mcParticles.iteratorAt(i);
      int mpdg = mother.pdgCode();
      if (abs(mpdg) == pdg && mpdg * p.pdgCode() > 0) { // check for quark
        if (quark_id > -1 || next_mother_id > -1) {     // we already found a possible candidate in the list of mothers, so now we have (at least) two
          // LOG(warning) << "Flavour tracking is ambiguous. Stopping here.";
          return -1;
        }
        quark_id = i;
      } else if ((static_cast<int>(abs(mpdg) / 100) == pdg || static_cast<int>(abs(mpdg) / 1000) == pdg) && mpdg * p.pdgCode() > 0) { // check for other mothers with flavour content
        if (quark_id > -1 || next_mother_id > -1) {                                                                                   // we already found a possible candidate in the list of mothers, so now we have (at least) two
          // LOG(warning) << "Flavour tracking is ambiguous. Stopping here.";
          return -1;
        }
        next_mother_id = i;
      }
    }
    if (quark_id > -1) { // we found the quark
      return quark_id;
    }
    if (next_mother_id > -1) {
      auto mother = mcParticles.iteratorAt(next_mother_id);
      return searchMothers(mother, mcParticles, pdg, equal); // we found at least something with flavour content, go deeper recursively
    }
    return -1;
  } else { // searching for first ancestor that is not quark anymore
    int quark_id = -1;
    for (int i : allmothersids) {
      auto mother = mcParticles.iteratorAt(i);
      int mpdg = abs(mother.pdgCode());
      if (mpdg == pdg && mother.pdgCode() == p.pdgCode()) { // found the quark
        if (quark_id > -1) {                                // we already found a possible candidate in the list of mothers, so now we have (at least) two
          // LOG(warning) << "Flavour tracking is ambiguous. Stopping here.";
          return -1;
        }
        quark_id = i;
      }
    }
    if (quark_id > -1) { // we found a quark and now search recursivly through entire history
      auto mother = mcParticles.iteratorAt(quark_id);
      return searchMothers(mother, mcParticles, pdg, equal);
    } else { // no HF quark as mother, so this is the ancestor
      return allmothersids[0];
    }
  }
  return -1;
}

template <typename T, typename U>
int findHFOrigin(T& p, U& mcParticles, int pdg)
{
  int quark_id = searchMothers(p, mcParticles, pdg, true); // try to find the hf quark
  if (quark_id == -1) {
    return -1;
  }
  auto quark = mcParticles.iteratorAt(quark_id);
  int id = searchMothers(quark, mcParticles, pdg, false); // try to find the first ancestor that is not the hf quark anymore
  return id;
}

template <typename T, typename U>
bool checkFromSameQuarkPair(T& p1, T& p2, U& mcParticles, int pdg)
{ // check if two particles come from the same hf q-qbar pair
  int id1 = findHFOrigin(p1, mcParticles, pdg);
  int id2 = findHFOrigin(p2, mcParticles, pdg);
  return id1 == id2 && id1 > -1 && id2 > -1;
}
//_______________________________________________________________________
template <typename T>
bool isCharmMeson(T const& track)
{
  if (400 < abs(track.pdgCode()) && abs(track.pdgCode()) < 500) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isCharmBaryon(T const& track)
{
  if (4000 < abs(track.pdgCode()) && abs(track.pdgCode()) < 5000) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isBeautyMeson(T const& track)
{
  if (500 < abs(track.pdgCode()) && abs(track.pdgCode()) < 600) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isBeautyBaryon(T const& track)
{
  if (5000 < abs(track.pdgCode()) && abs(track.pdgCode()) < 6000) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::utils::mcutil
#endif // PWGEM_DILEPTON_UTILS_MCUTILITIES_H_
