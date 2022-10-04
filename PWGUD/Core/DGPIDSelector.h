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

using namespace o2;

using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTrackCollisionIDs, aod::UDTracksPID, aod::UDTracksExtra>;
using UDTrackFull = UDTracksFull::iterator;

const int numDGPIDCutParameters = 9;
float particleMass(TDatabasePDG* pdg, int pid);

// -----------------------------------------------------------------------------
//  numPart:    Particle number to which the parameters apply
//  cutPID:     DPG code of particle hypothesis for nSigma calculation
//  cutDetector:detector: 1: TPC
//                        2: TOF
//  cutType:    cut type:  1: pt and nSigma within limits
//                        -1: nSigma out of limits within pt range
//                         2: pt and detector signal  within limits
//                        -2: detector signal out of limits within pt range
//  cutApply:   How to apply cut: 0: not active
//                                1: if information available
//                                2: return false if information not available
//  ptMin, ptMax, nSigmaTPCmin, nSigmaTPCmax, nSigmaTOFmin, nSigmaTOFmax: cut limits
struct DGPIDCut {

 public:
  // constructor
  DGPIDCut();
  DGPIDCut(float numPart, float cutPID, float cutDetector, float cutType, float cutApply,
           float ptMin, float ptMax, float nSigmamin, float nSigmamax);
  DGPIDCut(float* cutValues);
  ~DGPIDCut();

  // setters

  // getters
  void Print();
  int nPart() { return mnumPart; }
  int cutPID() { return mcutPID; }
  int cutDetector() { return mcutDetector; }
  int cutType() { return mcutType; }
  int cutApply() { return mcutApply; }
  float cutPtMin() { return mptMin; }
  float cutPtMax() { return mptMax; }
  float cutdetValueMin() { return mdetValuemin; }
  float cutdetValueMax() { return mdetValuemax; }

 private:
  int mnumPart;
  int mcutPID;
  int mcutDetector;
  int mcutType;
  int mcutApply;
  float mptMin;
  float mptMax;
  float mdetValuemin;
  float mdetValuemax;

  // ClassDefNV(DGPIDCut, 1);
};

// =============================================================================
// Cuts in mPIDCuts are combined (&&)
struct DGPIDCuts {
 public:
  // constructor
  DGPIDCuts();
  DGPIDCuts(std::vector<float> PIDCutValues);
  ~DGPIDCuts();

  // setter
  void clear()
  {
    mDGPIDCuts.clear();
  }
  void setPIDCuts(std::vector<float> PIDCutValues);

  // getter
  void Print();
  std::vector<DGPIDCut> Cuts() { return mDGPIDCuts; }

 private:
  std::vector<DGPIDCut> mDGPIDCuts;

  // ClassDefNV(DGPIDCuts, 1);
};

// =============================================================================
// object to hold customizable analysis parameters
struct DGAnaparHolder {
 public:
  // constructor
  DGAnaparHolder();
  DGAnaparHolder(int nCombine, std::vector<float> DGPIDs, std::vector<float> DGPIDCutValues);
  ~DGAnaparHolder();

  // getter
  void Print();
  int nCombine() const { return mNCombine; }
  std::vector<float> PIDs() { return mDGPIDs; }
  DGPIDCuts PIDCuts();
  std::vector<int> uniquePermutations();

 private:
  // helper functions
  void permutations(std::vector<uint>& ref, int n0, int np, std::vector<std::vector<uint>>& perms);
  int permutations(int n0, std::vector<std::vector<uint>>& perms);
  void makeUniquePermutations();

  // number of tracks to combine
  int mNCombine;

  // PID information
  std::vector<float> mDGPIDs;
  std::vector<float> mDGPIDCutValues;

  // unique permutations
  std::vector<int> muniquePerms;

  // ClassDefNV(DGAnaparHolder, 1);
};

// =============================================================================
// a structure which holds the indices of tracks and their invariant mass
struct DGParticle {
 public:
  DGParticle();
  DGParticle(TDatabasePDG* pdg, DGAnaparHolder anaPars, UDTracksFull const& tracks, std::vector<uint> comb);
  ~DGParticle();

  // getter
  void Print();
  std::vector<uint> trkinds() { return mtrkinds; }
  float M() { return mIVM.M(); }
  float Perp() { return mIVM.Perp(); }

 private:
  // invariant mass
  TLorentzVector mIVM;

  // indices of tracks included
  std::vector<uint> mtrkinds;

  // ClassDefNV(DGParticle, 1);
};

// =============================================================================
// A class to check PIDs of tracks and provide track combinations
struct DGPIDSelector {
 public:
  DGPIDSelector();
  ~DGPIDSelector();

  // setters
  void init(DGAnaparHolder anaPars);

  // getters
  void Print();
  bool isGoodTrack(UDTrackFull track, int cnt);
  int computeIVMs(UDTracksFull const& tracks);

  DGAnaparHolder getAnaPars() { return mAnaPars; }
  float getTPCnSigma(UDTrackFull track, int pid);
  float getTOFnSigma(UDTrackFull track, int pid);
  std::vector<DGParticle> IVMs() { return mIVMs; }

  int pid2ind(int pid);

 private:
  // analysis parameters
  DGAnaparHolder mAnaPars;

  // list of DGParticles
  std::vector<DGParticle> mIVMs;

  // particle properties
  TDatabasePDG* fPDG;

  // helper functions for computeIVMs
  void combinations(int n0, std::vector<uint>& pool, int np, std::vector<uint>& inds, int n,
                    std::vector<std::vector<uint>>& combs);
  int combinations(int n0, int np, std::vector<std::vector<uint>>& combs);
  std::vector<std::vector<uint>> combinations(int nPool);

  // ClassDefNV(DGPIDSelector, 1);
};

// -----------------------------------------------------------------------------
#endif // O2_ANALYSIS_DGPID_SELECTOR_
