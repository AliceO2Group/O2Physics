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

#ifndef O2_ANALYSIS_DGPID_SELECTOR_
#define O2_ANALYSIS_DGPID_SELECTOR_

#include <gandiva/projector.h>
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "PWGUD/DataModel/UDTables.h"
#include "DGAnaparHolder.h"

using namespace o2;

using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTrackCollisionIDs, aod::UDTracksPID>;
using UDTrackFull = UDTracksFull::iterator;

float particleMass(TDatabasePDG* pdg, int pid);

// -----------------------------------------------------------------------------
// a structure which holds the indices of tracks and their invariant mass
struct DGParticle {
 public:
  DGParticle() = default;
  DGParticle(TDatabasePDG* pdg, DGAnaparHolder anaPars, UDTracksFull const& tracks, std::vector<uint> comb);

  // getter
  std::vector<uint> trkinds() { return mtrkinds; }
  float M() { return mIVM.M(); }
  float Perp() { return mIVM.Perp(); }
  void Print();

 private:
  // invariant mass
  TLorentzVector mIVM;

  // indices of tracks included
  std::vector<uint> mtrkinds;

  ClassDefNV(DGParticle, 1);
};

// -----------------------------------------------------------------------------
// A class to check PIDs of tracks and provide track combinations
struct DGPIDSelector {
 public:
  DGPIDSelector();

  // setters
  void init(DGAnaparHolder anaPars)
  {
    mAnaPars = anaPars;
    mIVMs.clear();
  };

  // getters
  std::vector<DGParticle> IVMs() { return mIVMs; }
  float getTPCnSigma(UDTrackFull track, int hypo);
  bool isGoodTrack(UDTrackFull track, int cnt);
  int computeIVMs(int nCombine, UDTracksFull const& tracks);

 private:
  // analysis parameters
  DGAnaparHolder mAnaPars;

  // list of DGParticles
  std::vector<DGParticle> mIVMs;

  // particle properties
  TDatabasePDG* fPDG;
  int pid2ind(int pid);

  // helper functions for computeIVMs
  void permutations(std::vector<uint>& ref, int n0, int np, std::vector<std::vector<uint>>& perms);
  int permutations(int n0, std::vector<std::vector<uint>>& perms);
  void combinations(int n0, std::vector<uint>& pool, int np, std::vector<uint>& inds, int n,
                    std::vector<std::vector<uint>>& combs);
  int combinations(int n0, int np, std::vector<std::vector<uint>>& combs);
  std::vector<std::vector<uint>> combinations(int nCombine, int nPool);

  ClassDefNV(DGPIDSelector, 1);
};

// -----------------------------------------------------------------------------
#endif // O2_ANALYSIS_DGPID_SELECTOR_
