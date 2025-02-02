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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
/* Schematic Monte Carlo signal definition:

  Prong_0:      Particle(0,0) <-- Particle(0,1) <-- ... Prong(0,m_0) ... <-- Particle(0, n_0)
                                                          |
  Prong_1:      Particle(1,0) <-- Particle(1,1) <-- ... Prong(1,m_1) ... <-- Particle(1, n_1)
                                                          |
  Prong_i:      Particle(i,0) <-- Particle(i,1) <-- ... Prong(i,m_i) ... <-- Particle(i, n_i)
                                                          |
  Prong_N:      Particle(N,0) <-- Particle(N,1) <-- ... Prong(N,m_N) ... <-- Particle(N, n_N)

The MC signal model is composed of N prongs (see MCProng.h) with a variable number of generations specified for each prong (n_i, i=1,N).
At least one prong has to be specified for a valid MCsignal. An optional common ancestor for all prongs can be specified (m_i, i=1,N).
The common ancestor index, m_i, does not need to be the same for all prongs (e.g. some prongs can have a smaller number of generations in the history up to the common ancestor).
When a tuple of MC particles is checked whether it matches the specified signal, all prongs must be matched exactly with the particles in the tuple.
Example usage:

ATask {

MCSignal* mySignal;
MCSignal* mySignal2;

init() {
  MCProng prongElectronNonPromptJpsi(3,{11,443,502},{true,true,true},{false,false,false},{0,0,0},{0,0,0},{false,false,false});
  mySignal = new MCSignal("jpsiBeautyElectron", "Electrons from beauty jpsi decays", {prongElectronNonPromptJpsi}, {-1});
  MCProng prongAllFromBeauty(2,{0,503},{true,true},{false,false},{0,0},{0,0},{false,false});
  mySignal2 = new MCSignal("everythingFromBeautyPairs", "Everything from beauty hadrons pair", {prongAllFromBeauty,prongAllFromBeauty}, {-1,-1});
}
process(aod::McParticles const& mcTracks) {
  ...
  for (auto& mctrack : mcTracks) {
    if(mySignal.CheckSignal(true,mcTracks,mctrack)) {
      cout << "Found signal" << endl;
    }
  }
  for (auto& [mt1, mt2] : combinations(mcTracks, mcTracks)) {
    if(mySignal2.CheckSignal(true,mcTracks,mt1,mt2)) {
      cout << "Match found for " << mySignal2.GetName() << endl;
    }
  }
}
}
*/
#ifndef PWGDQ_CORE_MCSIGNAL_H_
#define PWGDQ_CORE_MCSIGNAL_H_

#include "MCProng.h"
#include "TNamed.h"

#include <vector>
#include <iostream>

class MCSignal : public TNamed
{
 public:
  MCSignal();
  MCSignal(int nProngs, const char* name = "", const char* title = ""); // NOLINT
  MCSignal(const char* name, const char* title, std::vector<MCProng> prongs, std::vector<int8_t> commonAncestors, bool excludeCommonAncestor = false);
  MCSignal(const MCSignal& c) = default;
  ~MCSignal() override = default;

  void SetProngs(std::vector<MCProng> prongs, std::vector<int8_t> commonAncestors);
  void AddProng(MCProng prong, int8_t commonAncestor = -1);
  void SetDecayChannelIsExclusive(int nProngs, bool option = true)
  {
    fDecayChannelIsExclusive = option;
    fNAncestorDirectProngs = nProngs;
  }
  void SetDecayChannelIsNotExclusive(int nProngs, bool option = true)
  {
    fDecayChannelIsNotExclusive = option;
    fNAncestorDirectProngs = nProngs;
  }

  int GetNProngs() const
  {
    return fNProngs;
  }
  int GetNGenerations() const
  {
    return fProngs[0].fNGenerations;
  }
  bool GetDecayChannelIsExclusive() const
  {
    return fDecayChannelIsExclusive;
  }
  bool GetDecayChannelIsNotExclusive() const
  {
    return fDecayChannelIsNotExclusive;
  }
  int GetNAncestorDirectProngs() const
  {
    return fNAncestorDirectProngs;
  }

  template <typename... T>
  bool CheckSignal(bool checkSources, const T&... args)
  {
    // Make sure number of tracks provided is equal to the number of prongs
    if (sizeof...(args) != fNProngs) { // TODO: addd a proper error message
      return false;
    }

    return CheckMC(0, checkSources, args...);
  };

  void PrintConfig();

 private:
  std::vector<MCProng> fProngs;            // vector of MCProng
  unsigned int fNProngs;                   // number of prongs
  std::vector<int8_t> fCommonAncestorIdxs; // index of the most recent ancestor, relative to each prong's history
  bool fExcludeCommonAncestor;             // explicitly request that there is no common ancestor
  bool fDecayChannelIsExclusive;           // if true, then the indicated mother particle has a number of daughters which is equal to the number of direct prongs defined in this MC signal
  bool fDecayChannelIsNotExclusive;        // if true, then the indicated mother particle has a number of daughters which is larger than the number of direct prongs defined in this MC signal
  int fNAncestorDirectProngs;              // number of direct prongs belonging to the common ancestor specified by this signal
  int fTempAncestorLabel;

  template <typename T>
  bool CheckProng(int i, bool checkSources, const T& track);

  bool CheckMC(int, bool)
  {
    return true;
  };

  template <typename T, typename... Ts>
  bool CheckMC(int i, bool checkSources, const T& track, const Ts&... args)
  {
    // recursive call of CheckMC for all args
    if (!CheckProng(i, checkSources, track)) {
      return false;
    } else {
      return CheckMC(i + 1, checkSources, args...);
    }
  };
};

template <typename T>
bool MCSignal::CheckProng(int i, bool checkSources, const T& track)
{
  using P = typename T::parent_t;
  auto currentMCParticle = track;

  // loop over the generations specified for this prong
  for (int j = 0; j < fProngs[i].fNGenerations; j++) {
    // check the PDG code
    if (!fProngs[i].TestPDG(j, currentMCParticle.pdgCode())) {
      return false;
    }
    // check the common ancestor (if specified)
    if (fNProngs > 1 && fCommonAncestorIdxs[i] == j) {
      if (i == 0) {
        fTempAncestorLabel = currentMCParticle.globalIndex();
        // In the case of decay channels marked as being "exclusive", check how many decay daughters this mother has registered
        //   in the stack and compare to the number of prongs defined for this MCSignal.
        //  If these numbers are equal, it means this decay MCSignal match is exclusive (there are no additional prongs for this mother besides the
        //     prongs defined here).
        if (currentMCParticle.has_daughters()) {
          if (fDecayChannelIsExclusive && currentMCParticle.daughtersIds()[1] - currentMCParticle.daughtersIds()[0] + 1 != fNAncestorDirectProngs) {
            return false;
          }
          if (fDecayChannelIsNotExclusive && currentMCParticle.daughtersIds()[1] - currentMCParticle.daughtersIds()[0] + 1 == fNAncestorDirectProngs) {
            return false;
          }
        }
      } else {
        if (currentMCParticle.globalIndex() != fTempAncestorLabel && !fExcludeCommonAncestor)
          return false;
        else if (currentMCParticle.globalIndex() == fTempAncestorLabel && fExcludeCommonAncestor)
          return false;
      }
    }

    // Update the currentMCParticle by moving either back in time (towards mothers, grandmothers, etc)
    // or in time (towards daughters) depending on how this was configured in the MC Signal
    if (!fProngs[i].fCheckGenerationsInTime) {
      // make sure that a mother exists in the stack before moving one generation further in history
      if (!currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
        return false;
      }
      if (currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
        currentMCParticle = currentMCParticle.template mothers_first_as<P>();
      }
    } else {
      // make sure that a daughter exists in the stack before moving one generation younger
      if (!currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
        return false;
      }
      if (currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
        const auto& daughtersSlice = currentMCParticle.template daughters_as<P>();
        for (auto& d : daughtersSlice) {
          if (fProngs[i].TestPDG(j + 1, d.pdgCode())) {
            currentMCParticle = d;
            break;
          }
        }
      }
    }
  }

  // check the various specified sources
  if (checkSources) {
    currentMCParticle = track;
    for (int j = 0; j < fProngs[i].fNGenerations; j++) {
      // check whether sources are required for this generation
      if (!fProngs[i].fSourceBits[j]) {
        continue;
      }
      // check each source
      uint64_t sourcesDecision = 0;
      // Check kPhysicalPrimary
      if (fProngs[i].fSourceBits[j] & (static_cast<uint64_t>(1) << MCProng::kPhysicalPrimary)) {
        if ((fProngs[i].fExcludeSource[j] & (static_cast<uint64_t>(1) << MCProng::kPhysicalPrimary)) != currentMCParticle.isPhysicalPrimary()) {
          sourcesDecision |= (static_cast<uint64_t>(1) << MCProng::kPhysicalPrimary);
        }
      }
      // Check kProducedInTransport
      if (fProngs[i].fSourceBits[j] & (static_cast<uint64_t>(1) << MCProng::kProducedInTransport)) {
        if ((fProngs[i].fExcludeSource[j] & (static_cast<uint64_t>(1) << MCProng::kProducedInTransport)) != (!currentMCParticle.producedByGenerator())) {
          sourcesDecision |= (static_cast<uint64_t>(1) << MCProng::kProducedInTransport);
        }
      }
      // Check kProducedByGenerator
      if (fProngs[i].fSourceBits[j] & (static_cast<uint64_t>(1) << MCProng::kProducedByGenerator)) {
        if ((fProngs[i].fExcludeSource[j] & (static_cast<uint64_t>(1) << MCProng::kProducedByGenerator)) != currentMCParticle.producedByGenerator()) {
          sourcesDecision |= (static_cast<uint64_t>(1) << MCProng::kProducedByGenerator);
        }
      }
      // Check kFromBackgroundEvent
      if (fProngs[i].fSourceBits[j] & (static_cast<uint64_t>(1) << MCProng::kFromBackgroundEvent)) {
        if ((fProngs[i].fExcludeSource[j] & (static_cast<uint64_t>(1) << MCProng::kFromBackgroundEvent)) != currentMCParticle.fromBackgroundEvent()) {
          sourcesDecision |= (static_cast<uint64_t>(1) << MCProng::kFromBackgroundEvent);
        }
      }
      // no source bit is fulfilled
      if (!sourcesDecision) {
        return false;
      }
      // if fUseANDonSourceBitMap is on, request all bits
      if (fProngs[i].fUseANDonSourceBitMap[j] && (sourcesDecision != fProngs[i].fSourceBits[j])) {
        return false;
      }

      // Update the currentMCParticle by moving either back in time (towards mothers, grandmothers, etc)
      //    or in time (towards daughters) depending on how this was configured in the MSignal
      if (!fProngs[i].fCheckGenerationsInTime) {
        // make sure that a mother exists in the stack before moving one generation further in history
        if (!currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
          return false;
        }
        if (currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
          currentMCParticle = currentMCParticle.template mothers_first_as<P>();
        }
      } else {
        // prong history will be moved to the branch of the first daughter that matches the PDG requirement
        // make sure that a daughter exists in the stack before moving one generation younger
        if (!currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
          return false;
        }
        if (currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
          const auto& daughtersSlice = currentMCParticle.template daughters_as<P>();
          for (auto& d : daughtersSlice) {
            if (fProngs[i].TestPDG(j + 1, d.pdgCode())) {
              currentMCParticle = d;
              break;
            }
          }
        }
      }
    }
  }

  if (fProngs[i].fPDGInHistory.size() == 0) {
    return true;
  } else { // check if mother pdg is in history
    std::vector<int> pdgInHistory;

    // while find mothers, check if the provided PDG codes are included or excluded in the particle decay history
    unsigned int nIncludedPDG = 0;
    for (unsigned int k = 0; k < fProngs[i].fPDGInHistory.size(); k++) {
      currentMCParticle = track;
      if (!fProngs[i].fExcludePDGInHistory[k])
        nIncludedPDG++;
      int ith = 0;

      // Note: Currently no need to check generation InTime, so disable if case and always check BackInTime (direction of mothers)
      //       The option to check for daughter in decay chain is still implemented but commented out.

      // if (!fProngs[i].fCheckGenerationsInTime) { // check generation back in time
      while (currentMCParticle.has_mothers()) {
        auto mother = currentMCParticle.template mothers_first_as<P>();
        if (!fProngs[i].fExcludePDGInHistory[k] && fProngs[i].ComparePDG(mother.pdgCode(), fProngs[i].fPDGInHistory[k], true, fProngs[i].fExcludePDGInHistory[k])) {
          pdgInHistory.emplace_back(mother.pdgCode());
          break;
        }
        if (fProngs[i].fExcludePDGInHistory[k] && !fProngs[i].ComparePDG(mother.pdgCode(), fProngs[i].fPDGInHistory[k], true, fProngs[i].fExcludePDGInHistory[k])) {
          return false;
        }
        ith++;
        currentMCParticle = mother;
        if (ith > 10) { // need error message. Given pdg code was not found within 10 generations of the particles decay chain.
          break;
        }
      }
      // } else { // check generation in time
      //   if (!currentMCParticle.has_daughters())
      //     return false;
      //   const auto& daughtersSlice = currentMCParticle.template daughters_as<P>();
      //   for (auto& d : daughtersSlice) {
      //     if (!fProngs[i].fExcludePDGInHistory[k] && fProngs[i].ComparePDG(d.pdgCode(), fProngs[i].fPDGInHistory[k], true, fProngs[i].fExcludePDGInHistory[k])) {
      //       pdgInHistory.emplace_back(d.pdgCode());
      //       break;
      //     }
      //     if (fProngs[i].fExcludePDGInHistory[k] && !fProngs[i].ComparePDG(d.pdgCode(), fProngs[i].fPDGInHistory[k], true, fProngs[i].fExcludePDGInHistory[k])) {
      //       return false;
      //     }
      //     ith++;
      //     if (ith > 10) { // need error message. Given pdg code was not found within 10 generations of the particles decay chain.
      //       break;
      //     }
      //   }
      // }
    }
    if (pdgInHistory.size() != nIncludedPDG) { // vector has as many entries as mothers (daughters) defined for prong
      return false;
    }
  }
  return true;
}

#endif // PWGDQ_CORE_MCSIGNAL_H_
