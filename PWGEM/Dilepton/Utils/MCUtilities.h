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

#include "Framework/AnalysisDataModel.h"
#include "Framework/Logger.h"

#include <algorithm>
#include <string>
#include <vector>

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
template <typename TTrack>
int hasFakeMatchITSTPC(TTrack const& track)
{
  // track and mctracklabel have to be joined.
  // bit 13 -- ITS/TPC labels are not equal

  if ((track.mcMask() & 1 << 13)) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename TTrack>
int hasFakeMatchITSTPCTOF(TTrack const&)
{
  // track and mctracklabel have to be joined.
  return false;
  // if ((track.mcMask() & 1 << 13) && (track.mcMask() & 1 << 15)) {
  //   return true;
  // } else {
  //   return false;
  // }
}
//_______________________________________________________________________
template <typename TTrack>
int hasFakeMatchMFTMCH(TTrack const& track)
{
  // fwdtrack and mcfwdtracklabel have to be joined.
  if ((track.mcMask() & 1 << 7)) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isCharmonia(T const& track)
{
  if (std::abs(track.pdgCode()) < 100) {
    return false;
  }

  std::string pdgStr = std::to_string(std::abs(track.pdgCode()));
  int n = pdgStr.length();
  int pdg3 = std::stoi(pdgStr.substr(n - 3, 3));

  if (pdg3 == 441 || pdg3 == 443 || pdg3 == 445 || pdg3 == 447) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isCharmMeson(T const& track)
{
  if (isCharmonia(track)) {
    return false;
  }

  if (400 < std::abs(track.pdgCode()) && std::abs(track.pdgCode()) < 500) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isCharmBaryon(T const& track)
{
  if (4000 < std::abs(track.pdgCode()) && std::abs(track.pdgCode()) < 5000) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isBottomonia(T const& track)
{
  if (std::abs(track.pdgCode()) < 100) {
    return false;
  }

  std::string pdgStr = std::to_string(std::abs(track.pdgCode()));
  int n = pdgStr.length();
  int pdg3 = std::stoi(pdgStr.substr(n - 3, 3));

  if (pdg3 == 551 || pdg3 == 553 || pdg3 == 555 || pdg3 == 557) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isBeautyMeson(T const& track)
{
  if (isBottomonia(track)) {
    return false;
  }

  if (500 < std::abs(track.pdgCode()) && std::abs(track.pdgCode()) < 600) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T>
bool isBeautyBaryon(T const& track)
{
  if (5000 < std::abs(track.pdgCode()) && std::abs(track.pdgCode()) < 6000) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T, typename U>
bool isWeakDecayFromCharmHadron(T const& mcParticle, U const& mcParticles)
{
  // require that the direct mother is charm hadron via semileptonic. e.g. hc->e, not hc->X->pi0->eegamma
  if (!mcParticle.has_mothers()) {
    return false;
  }
  // if (mcParticle.getProcess() != 4) { // weak decay
  //   return false;
  // }
  auto mp = mcParticles.iteratorAt(mcParticle.mothersIds()[0]);
  if (isCharmMeson(mp) || isCharmBaryon(mp)) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename T, typename U>
bool isWeakDecayFromBeautyHadron(T const& mcParticle, U const& mcParticles)
{
  // require that the direct mother is beauty hadron via semileptonice decay. e.g. hb->e, not hb->X->pi0->eegamma
  if (!mcParticle.has_mothers()) {
    return false;
  }
  // if (mcParticle.getProcess() != 4) { // weak decay
  //   return false;
  // }
  auto mp = mcParticles.iteratorAt(mcParticle.mothersIds()[0]);
  if (isBeautyMeson(mp) || isBeautyBaryon(mp)) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
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
int getMotherPDGCode(TMCParticle const& p, TMCParticles const& mcparticles)
{
  int motherid = p.mothersIds()[0];
  auto mother = mcparticles.iteratorAt(motherid);
  return (mother.pdgCode());
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
  if (std::abs(mp_tmp.pdgCode()) < 1e+9 && (std::to_string(std::abs(mp_tmp.pdgCode()))[std::to_string(std::abs(mp_tmp.pdgCode())).length() - 2] == '5' && std::to_string(std::abs(mp_tmp.pdgCode()))[std::to_string(std::abs(mp_tmp.pdgCode())).length() - 3] == '5') && std::abs(mp_tmp.pdgCode()) % 2 == 1) {
    return -999; // reject bottomonia
  }

  while (motherid > -1) {
    if (motherid < mcparticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcparticles.iteratorAt(motherid);
      if (std::abs(mp.pdgCode()) < 1e+9 && (std::to_string(std::abs(mp.pdgCode()))[std::to_string(std::abs(mp.pdgCode())).length() - 3] == '5' || std::to_string(std::abs(mp.pdgCode()))[std::to_string(std::abs(mp.pdgCode())).length() - 4] == '5')) {
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

//_______________________________________________________________________
template <typename TMCParticle, typename TMCParticles>
int IsFromCharm(TMCParticle const& p, TMCParticles const& mcparticles)
{
  if (!p.has_mothers()) {
    return -999;
  }

  int motherid = p.mothersIds()[0]; // first mother index
  auto mp_tmp = mcparticles.iteratorAt(motherid);
  if (std::abs(mp_tmp.pdgCode()) < 1e+9 && (std::to_string(std::abs(mp_tmp.pdgCode()))[std::to_string(std::abs(mp_tmp.pdgCode())).length() - 2] == '4' && std::to_string(std::abs(mp_tmp.pdgCode()))[std::to_string(std::abs(mp_tmp.pdgCode())).length() - 3] == '4') && std::abs(mp_tmp.pdgCode()) % 2 == 1) {
    return -999; // reject bottomonia
  }
  while (motherid > -1) {
    if (motherid < mcparticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcparticles.iteratorAt(motherid);
      if (std::abs(mp.pdgCode()) < 1e+9 && (std::to_string(std::abs(mp.pdgCode()))[std::to_string(std::abs(mp.pdgCode())).length() - 3] == '4' || std::to_string(std::abs(mp.pdgCode()))[std::to_string(std::abs(mp.pdgCode())).length() - 4] == '4')) {
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
//_______________________________________________________________________
template <typename TMCParticle, typename TMCParticles>
bool findFlavorOscillationB(TMCParticle const& mcParticle, TMCParticles const& mcParticles)
{
  // B0 or B0s can oscillate.
  if (!mcParticle.has_mothers()) {
    return false;
  }

  // store all mother1 relation
  int motherid1 = mcParticle.mothersIds()[0]; // first mother index
  while (motherid1 > -1) {
    if (motherid1 < mcParticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcParticles.iteratorAt(motherid1);
      if ((std::abs(mp.pdgCode()) == 511 || std::abs(mp.pdgCode()) == 531) && mp.getGenStatusCode() == 92) {
        return true;
      }

      if (mp.has_mothers()) {
        motherid1 = mp.mothersIds()[0];
      } else {
        motherid1 = -999;
      }
    } else {
      LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid1, mcParticles.size());
    }
  }
  return false;
}
//_______________________________________________________________________
template <typename TMCParticle, typename TMCParticles>
int find1stHadron(TMCParticle const& mcParticle, TMCParticles const& mcParticles)
{
  // find 1st hadron in decay chain except beam.
  if (!mcParticle.has_mothers()) {
    return -1;
  }

  // store all mother1 relation
  std::vector<int> mothers_id;
  std::vector<int> mothers_pdg;

  int motherid1 = mcParticle.mothersIds()[0]; // first mother index
  while (motherid1 > -1) {
    if (motherid1 < mcParticles.size()) { // protect against bad mother indices. why is this needed?
      auto mp = mcParticles.iteratorAt(motherid1);
      mothers_id.emplace_back(motherid1);
      mothers_pdg.emplace_back(mp.pdgCode());

      if (mp.has_mothers()) {
        motherid1 = mp.mothersIds()[0];
      } else {
        motherid1 = -999;
      }
    } else {
      LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid1, mcParticles.size());
    }
  }

  int counter = 0;
  for (const auto& pdg : mothers_pdg) {
    if (std::abs(pdg) <= 6 || std::abs(pdg) == 21 || (std::abs(pdg) == 2212 && counter == static_cast<int>(mothers_pdg.size() - 1)) || (std::abs(pdg) > 1e+9 && counter == static_cast<int>(mothers_pdg.size() - 1))) { // quarks or gluon or  proton or ion beam
      break;
    }
    counter++;
  }

  int hadronId = -1;
  if (counter == 0) { // particle directly from beam // only for protection.
    hadronId = mcParticle.globalIndex();
  } else {
    hadronId = mothers_id[counter - 1];
  }

  mothers_id.clear();
  mothers_id.shrink_to_fit();
  mothers_pdg.clear();
  mothers_pdg.shrink_to_fit();
  return hadronId;
}
//_______________________________________________________________________
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
      // LOGF(info, "mp1.globalIndex() = %d, mp1.pdgCode() = %d, mp1.getGenStatusCode() = %d", mp.globalIndex(), mp.pdgCode(), mp.getGenStatusCode());

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
      // LOGF(info, "mp2.globalIndex() = %d, mp2.pdgCode() = %d, mp2.getGenStatusCode() = %d", mp.globalIndex(), mp.pdgCode(), mp.getGenStatusCode());

      if (mp.has_mothers()) {
        motherid2 = mp.mothersIds()[0];
      } else {
        motherid2 = -999;
      }
    } else {
      LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid2, mcparticles.size());
    }
  }

  // require correlation between q-qbar. (not q-q) // need statusCode

  auto mpfh1 = mcparticles.iteratorAt(find1stHadron(p1, mcparticles));
  auto mpfh2 = mcparticles.iteratorAt(find1stHadron(p2, mcparticles));
  bool isFOFound1 = findFlavorOscillationB(p1, mcparticles);
  bool isFOFound2 = findFlavorOscillationB(p2, mcparticles);

  bool is_direct_from_b1 = isWeakDecayFromBeautyHadron(p1, mcparticles);
  bool is_direct_from_b2 = isWeakDecayFromBeautyHadron(p2, mcparticles);
  bool is_prompt_c1 = isWeakDecayFromCharmHadron(p1, mcparticles) && IsFromBeauty(p1, mcparticles) < 0;
  bool is_prompt_c2 = isWeakDecayFromCharmHadron(p2, mcparticles) && IsFromBeauty(p2, mcparticles) < 0;
  bool is_c_from_b1 = isWeakDecayFromCharmHadron(p1, mcparticles) && IsFromBeauty(p1, mcparticles) > 0;
  bool is_c_from_b2 = isWeakDecayFromCharmHadron(p2, mcparticles) && IsFromBeauty(p2, mcparticles) > 0;

  if (is_prompt_c1 && is_prompt_c2 && mpfh1.pdgCode() * mpfh2.pdgCode() < 0) { // charmed mesons never oscillate. only ULS
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

  bool b2l_b2l_case0 = is_direct_from_b1 && is_direct_from_b2 && mpfh1.pdgCode() * mpfh2.pdgCode() < 0 && !isFOFound1 && !isFOFound2;                 // 0 oscillation: bbbar -> ll ULS
  bool b2l_b2l_case1 = is_direct_from_b1 && is_direct_from_b2 && mpfh1.pdgCode() * mpfh2.pdgCode() > 0 && static_cast<bool>(isFOFound1 ^ isFOFound2); // 1 oscillation: bbbar -> ll LS
  bool b2l_b2l_case2 = is_direct_from_b1 && is_direct_from_b2 && mpfh1.pdgCode() * mpfh2.pdgCode() < 0 && isFOFound1 && isFOFound2;                   // 2 oscillation: bbbar -> ll ULS

  if (b2l_b2l_case0 || b2l_b2l_case1 || b2l_b2l_case2) {
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

  bool b2c2l_b2c2l_case0 = is_c_from_b1 && is_c_from_b2 && mpfh1.pdgCode() * mpfh2.pdgCode() < 0 && !isFOFound1 && !isFOFound2;                 // 0 oscillation: bbbar -> ll ULS
  bool b2c2l_b2c2l_case1 = is_c_from_b1 && is_c_from_b2 && mpfh1.pdgCode() * mpfh2.pdgCode() > 0 && static_cast<bool>(isFOFound1 ^ isFOFound2); // 1 oscillation: bbbar -> ll LS
  bool b2c2l_b2c2l_case2 = is_c_from_b1 && is_c_from_b2 && mpfh1.pdgCode() * mpfh2.pdgCode() < 0 && isFOFound1 && isFOFound2;                   // 2 oscillation: bbbar -> ll ULS

  if (b2c2l_b2c2l_case0 || b2c2l_b2c2l_case1 || b2c2l_b2c2l_case2) {
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
    // No pair sign oscillation due to B0(s) oscillation for the same mother.
    for (const auto& mid1 : mothers_id1) {
      for (const auto& mid2 : mothers_id2) {
        if (mid1 == mid2) {
          auto common_mp = mcparticles.iteratorAt(mid1);
          int mp_pdg = common_mp.pdgCode();
          bool is_mp_diquark = (1100 < std::abs(mp_pdg) && std::abs(mp_pdg) < 5600) && std::to_string(mp_pdg)[std::to_string(mp_pdg).length() - 2] == '0';
          if (!is_mp_diquark && std::abs(mp_pdg) < 1e+9 && (std::to_string(std::abs(mp_pdg))[std::to_string(std::abs(mp_pdg)).length() - 3] == '5' || std::to_string(std::abs(mp_pdg))[std::to_string(std::abs(mp_pdg)).length() - 4] == '5')) {
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
      } // end of motherid2
    } // end of motherid1

    bool b2c2l_b2l_diffb_case0 = mpfh1.pdgCode() * mpfh2.pdgCode() > 0 && !isFOFound1 && !isFOFound2;                 // 0 oscillation: bbbar -> ll LS
    bool b2c2l_b2l_diffb_case1 = mpfh1.pdgCode() * mpfh2.pdgCode() < 0 && static_cast<bool>(isFOFound1 ^ isFOFound2); // 1 oscillation: bbbar -> ll ULS
    bool b2c2l_b2l_diffb_case2 = mpfh1.pdgCode() * mpfh2.pdgCode() > 0 && isFOFound1 && isFOFound2;                   // 2 oscillation: bbbar -> ll LS

    if (b2c2l_b2l_diffb_case0 || b2c2l_b2l_diffb_case1 || b2c2l_b2l_diffb_case2) {
      mothers_id1.clear();
      mothers_pdg1.clear();
      mothers_id2.clear();
      mothers_pdg2.clear();
      mothers_id1.shrink_to_fit();
      mothers_pdg1.shrink_to_fit();
      mothers_id2.shrink_to_fit();
      mothers_pdg2.shrink_to_fit();
      return static_cast<int>(EM_HFeeType::kBCe_Be_DiffB); // b->c->e and b->e, decay type = 4. this should happen only in LS.
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
    } else if ((80 < std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) && std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) < 90) || (100 < std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) && std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) < 110)) { // NOTE: THIS IS GENERATOR DEPENDENT AND WORKS ONLY FOR PYTHIA!
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
      // if (std::abs(mpdg) == pdg && mpdg * p.pdgCode() > 0) { // check for quark
      if (std::abs(mpdg) == pdg) {                  // check for quark to allow for beauty and charm + oscillation
        if (quark_id > -1 || next_mother_id > -1) { // we already found a possible candidate in the list of mothers, so now we have (at least) two
          // LOG(warning) << "Flavour tracking is ambiguous. Stopping here.";
          return -1;
        }
        quark_id = i;
        //} else if ((static_cast<int>(std::abs(mpdg) / 100) == pdg || static_cast<int>(std::abs(mpdg) / 1000) == pdg) && mpdg * p.pdgCode() > 0) { // check for other mothers with flavour content
      } else if ((static_cast<int>(std::abs(mpdg) / 100) == pdg || static_cast<int>(std::abs(mpdg) / 1000) == pdg)) { // check for other mothers with flavour content to allow for beauty and charm
        if (quark_id > -1 || next_mother_id > -1) {                                                                   // we already found a possible candidate in the list of mothers, so now we have (at least) two
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
      int mpdg = std::abs(mother.pdgCode());
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

//_______________________________________________________________________
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
//_______________________________________________________________________
template <typename T, typename U>
bool checkFromSameQuarkPair(T& p1, T& p2, U& mcParticles, int pdg)
{ // check if two particles come from the same hf q-qbar pair
  int id1 = findHFOrigin(p1, mcParticles, pdg);
  int id2 = findHFOrigin(p2, mcParticles, pdg);
  return id1 == id2 && id1 > -1 && id2 > -1;
}
//_______________________________________________________________________
//_______________________________________________________________________
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::utils::mcutil
#endif // PWGEM_DILEPTON_UTILS_MCUTILITIES_H_
