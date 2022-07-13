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
#ifndef MCSignal_H
#define MCSignal_H

#include "MCProng.h"
#include "TNamed.h"

#include <vector>
#include <iostream>
using std::cout;
using std::endl;

class MCSignal : public TNamed
{
 public:
  MCSignal();
  MCSignal(int nProngs, const char* name = "", const char* title = "");
  MCSignal(const char* name, const char* title, std::vector<MCProng> prongs, std::vector<short> commonAncestors);
  MCSignal(const MCSignal& c) = default;
  ~MCSignal() override = default;

  void SetProngs(std::vector<MCProng> prongs, std::vector<short> commonAncestors);
  void AddProng(MCProng prong, short commonAncestor = -1);

  int GetNProngs() const
  {
    return fNProngs;
  }
  int GetNGenerations() const
  {
    return fProngs[0].fNGenerations;
  }

  template <typename U, typename... T>
  bool CheckSignal(bool checkSources, const U& mcStack, const T&... args)
  {
    // Make sure number of tracks provided is equal to the number of prongs
    if (sizeof...(args) != fNProngs) { // TODO: addd a proper error message
      return false;
    }

    return CheckMC(0, checkSources, mcStack, args...);
  };

  void PrintConfig();

 private:
  std::vector<MCProng> fProngs;
  unsigned int fNProngs;
  std::vector<short> fCommonAncestorIdxs;
  int fTempAncestorLabel;

  template <typename U, typename T>
  bool CheckProng(int i, bool checkSources, const U& mcStack, const T& track);

  template <typename U>
  bool CheckMC(int, bool, U)
  {
    return true;
  };

  template <typename U, typename T, typename... Ts>
  bool CheckMC(int i, bool checkSources, const U& mcStack, const T& track, const Ts&... args)
  {
    // recursive call of CheckMC for all args
    if (!CheckProng(i, checkSources, mcStack, track)) {
      return false;
    } else {
      return CheckMC(i + 1, checkSources, mcStack, args...);
    }
  };

  ClassDef(MCSignal, 1);
};

template <typename U, typename T>
bool MCSignal::CheckProng(int i, bool checkSources, const U& mcStack, const T& track)
{
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
      } else {
        if (currentMCParticle.globalIndex() != fTempAncestorLabel) {
          return false;
        }
      }
    }

    // if checking back in time: look for mother
    // else (checking in time): look for daughter
    if (!fProngs[i].fCheckGenerationsInTime) {
      // make sure that a mother exists in the stack before moving one generation further in history
      if (!currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
        return false;
      }
      /*for (auto& m : mcParticle.mothers_as<aod::McParticles_001>()) {
        LOGF(debug, "M2 %d %d", mcParticle.globalIndex(), m.globalIndex());
      }*/
      if (currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
        currentMCParticle = currentMCParticle.template mothers_first_as<U>();
        // currentMCParticle = mcStack.iteratorAt(currentMCParticle.mothersIds()[0]);
        // currentMCParticle = currentMCParticle.template mother0_as<U>();
      }
    } else {
      // make sure that a daughter exists in the stack before moving one generation younger
      if (!currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
        return false;
      }
      if (currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
        const auto& daughtersSlice = currentMCParticle.template daughters_as<U>();
        for (auto& d : daughtersSlice) {
          if (fProngs[i].TestPDG(j + 1, d.pdgCode())) {
            currentMCParticle = d;
            break;
          }
        }
      }
    }
  }

  if (checkSources) {
    currentMCParticle = track;
    for (int j = 0; j < fProngs[i].fNGenerations; j++) {
      if (!fProngs[i].fSourceBits[j]) {
        // no sources required for this generation
        continue;
      }
      // check each source
      uint64_t sourcesDecision = 0;
      // Check kPhysicalPrimary
      if (fProngs[i].fSourceBits[j] & (uint64_t(1) << MCProng::kPhysicalPrimary)) {
        if ((fProngs[i].fExcludeSource[j] & (uint64_t(1) << MCProng::kPhysicalPrimary)) != currentMCParticle.isPhysicalPrimary()) {
          sourcesDecision |= (uint64_t(1) << MCProng::kPhysicalPrimary);
        }
      }
      // Check kProducedInTransport
      if (fProngs[i].fSourceBits[j] & (uint64_t(1) << MCProng::kProducedInTransport)) {
        if ((fProngs[i].fExcludeSource[j] & (uint64_t(1) << MCProng::kProducedInTransport)) != (!currentMCParticle.producedByGenerator())) {
          sourcesDecision |= (uint64_t(1) << MCProng::kProducedInTransport);
        }
      }
      // Check kProducedByGenerator
      if (fProngs[i].fSourceBits[j] & (uint64_t(1) << MCProng::kProducedByGenerator)) {
        if ((fProngs[i].fExcludeSource[j] & (uint64_t(1) << MCProng::kProducedByGenerator)) != currentMCParticle.producedByGenerator()) {
          sourcesDecision |= (uint64_t(1) << MCProng::kProducedByGenerator);
        }
      }
      // Check kFromBackgroundEvent
      if (fProngs[i].fSourceBits[j] & (uint64_t(1) << MCProng::kFromBackgroundEvent)) {
        if ((fProngs[i].fExcludeSource[j] & (uint64_t(1) << MCProng::kFromBackgroundEvent)) != currentMCParticle.fromBackgroundEvent()) {
          sourcesDecision |= (uint64_t(1) << MCProng::kFromBackgroundEvent);
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

      // if checking back in time: look for mother
      // else (checking in time): look for daughter
      if (!fProngs[i].fCheckGenerationsInTime) {
        // move one generation back in history
        // make sure that a mother exists in the stack before moving one generation further in history
        if (!currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
          return false;
        }
        if (currentMCParticle.has_mothers() && j < fProngs[i].fNGenerations - 1) {
          /*for (auto& m : mcParticle.mothers_as<aod::McParticles_001>()) {
            LOGF(debug, "M2 %d %d", mcParticle.globalIndex(), m.globalIndex());
          }*/
          currentMCParticle = currentMCParticle.template mothers_first_as<U>();
          // currentMCParticle = mcStack.iteratorAt(currentMCParticle.mothersIds()[0]);
          // currentMCParticle = currentMCParticle.template mother0_as<U>();
        }
        /*if (j < fProngs[i].fNGenerations - 1) {
          currentMCParticle = mcStack.iteratorAt(currentMCParticle.mother0Id());
          //currentMCParticle = currentMCParticle.template mother0_as<U>();
        }*/
      } else {
        // move one generation further in history
        // pong history will be moved to the branch of the first daughter that matches the PDG requirement
        // make sure that a daughter exists in the stack before moving one generation younger
        if (!currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
          return false;
        }
        if (currentMCParticle.has_daughters() && j < fProngs[i].fNGenerations - 1) {
          const auto& daughtersSlice = currentMCParticle.template daughters_as<U>();
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

  return true;
}

#endif
