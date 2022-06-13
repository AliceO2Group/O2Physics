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

#ifndef O2_ANALYSIS_PID_SELECTOR_
#define O2_ANALYSIS_PID_SELECTOR_

#include "PWGUD/DataModel/DGCandidates.h"
#include "anaparHolder.h"

using namespace o2;

// -----------------------------------------------------------------------------
// a structure which holds the indices of tracks and their invariant mass
struct DGParticle {
 public:
  DGParticle() = default;
  DGParticle(anaparHolder anaPars, aod::DGTracks const& dgtracks, std::vector<uint> comb);

  // getter
  std::vector<uint> trkinds() { return mtrkinds; }
  float M() { return mM; }
  float Perp() { return mPerp; }
  void Print();

 private:
  // the invariant mass
  float mM;

  // pt of system
  float mPerp;

  // indices of tracks included
  std::vector<uint> mtrkinds;

  ClassDefNV(DGParticle, 1);
};

// -----------------------------------------------------------------------------
// A class to check PIDs of tracks and provide track combinations
struct pidSelector {
 public:
  pidSelector() = default;

  // setters
  void init(anaparHolder anaPars)
  {
    mAnaPars = anaPars;
    mIVMs.clear();
  };

  int computeIVMs(int nCombine, aod::DGTracks const& dgtracks);

  // getters
  std::vector<DGParticle> IVMs() { return mIVMs; }
  float getTPCnSigma(aod::DGTrack track, int hypo);
  bool isGoodTrack(aod::DGTrack track, int cnt);

 private:
  // analysis parameters
  anaparHolder mAnaPars;

  // list of DGParticles
  std::vector<DGParticle> mIVMs;

  ClassDefNV(pidSelector, 1);
};

// -----------------------------------------------------------------------------
#endif // O2_ANALYSIS_PID_SELECTOR_
